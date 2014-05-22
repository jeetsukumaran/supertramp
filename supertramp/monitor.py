#! /usr/bin/env python

import pandas

class AttributeTracker(object):

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
            record["prev_" + self.field_name] = getattr(simulator, "prev_" + self.attr_name)
            record["delta_" + self.field_name] = record[self.field_name] - record["prev_" + self.field_name]
            if simulator.current_gen == simulator.prev_sampled_generation:
                record["rate_of_change_" + self.field_name] = "NA"
            else:
                record["rate_of_change_" + self.field_name] = record["delta_" + self.field_name]  / (simulator.current_gen - simulator.prev_sampled_generation)
            setattr(simulator, "prev_" + self.attr_name, record[self.field_name])
        return record

class SimulatorMonitor(object):

    def __init__(self):
        self.trackers = []
        self.records = []

    def add_attribute_tracker(self, attr_name, field_name, sample_diffs):
        s = AttributeTracker(
                attr_name=attr_name,
                field_name=field_name,
                sample_diffs=sample_diffs,
                )
        self.trackers.append(s)
        return s

    def sample(self, simulator):
        record = {}
        if not hasattr(simulator, "prev_sampled_generation"):
            setattr(simulator, "prev_sampled_generation", 0)
        for tracker in self.trackers:
            tracker.sample(simulator, record)
        self.records.append(record)
        setattr(simulator, "prev_sampled_generation", simulator.current_gen)

    def as_data_frame(self):
        return pandas.DataFrame(self.records)


