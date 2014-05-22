#! /usr/bin/env python

import pandas

class AttributeSampler(object):

    def __init__(self,
            attr_name,
            field_name=None,
            sample_diffs=True):
        self.attr_name = attr_name
        if field_name is None:
            self.field_name = self.attr_name
        else:
            self.field_name = field_name
        self.sample_diffs = sample_diffs

    def sample(self, simulator, record):
        record[self.field_name] = float(getattr(simulator, self.attr_name))
        if self.sample_diffs:
            if not hasattr(simulator, "prev_" + self.attr_name):
                setattr(simulator, "prev_" + self.attr_name, 0.0)
            if not hasattr(simulator, "prev_sampled_generation"):
                setattr(simulator, "prev_sampled_generation", 0)
            record["prev_" + self.field_name] = "x"
            record["prev_" + self.field_name] = getattr(simulator, "prev_" + self.attr_name)
            record["delta_" + self.field_name] = record[self.field_name] - record["prev_" + self.field_name]
            if simulator.current_gen == simulator.prev_sampled_generation:
                record["rate_of_change_" + self.field_name] = "NA"
            else:
                record["rate_of_change_" + self.field_name] = record["delta_" + self.field_name]  / (simulator.current_gen - simulator.prev_sampled_generation)
        return record


class SampleRecord(dict):

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

    def sample_params(self, simulator):
        self.data["sampled_gen"] = simulator.current_gen
        self.data["s0"] = simulator.diversification_model_s0
        self.data["e0"] = simulator.diversification_model_e0

    def sample_data(self, simulator):
        # number of extant lineages
        # number of lineages produced
        # number of lineages extinct
        self.data["gen"] = simulator.current_gen
        for attr_name in (
                "num_extant_lineages",
                "num_births",
                "num_extinctions",
                "num_extirpations",
                ):
            self.sample_attribute(simulator, attr_name)

    def sample_attribute(self,
            simulator,
            attr_name,
            param_field_name=None):
        if param_field_name is None:
            param_field_name = attr_name
        self.data[param_field_name] = float(getattr(simulator, attr_name))
        if not hasattr(simulator, "prev_" + attr_name):
            setattr(simulator, "prev_" + attr_name, 0.0)
        if not hasattr(simulator, "prev_sampled_generation"):
            setattr(simulator, "prev_sampled_generation", 0)
        self.data["prev_" + param_field_name] = "x"
        self.data["prev_" + param_field_name] = getattr(simulator, "prev_" + attr_name)
        self.data["delta_" + param_field_name] = self.data[param_field_name] - self.data["prev_" + param_field_name]
        if simulator.current_gen == simulator.prev_sampled_generation:
            self.data["rate_of_change_" + param_filed_name] = "NA"
        else:
            self.data["rate_of_change_" + param_field_name] = self.data["delta_" + param_field_name]  / (simulator.current_gen - simulator.prev_sampled_generation)

    def as_dict(self):
        return dict(self)

class SimulatorSampler(object):

    def __init__(self):
        self.samplers = []
        self.samples = []

    def add_sampler(self, attr_name, field_name, sample_diffs):
        s = AttributeSampler(
                attr_name=attr_name,
                field_name=field_name,
                sample_diffs=sample_diffs,
                )
        self.samplers.append(s)
        return s

    def sample(self, simulator):
        record = SampleRecord()
        for sampler in self.samplers:
            sampler.sample(simulator, record)
        self.samples.append(record)

    def as_data_frame(self):
        return pandas.DataFrame([s.as_dict() for s in self.samples])


