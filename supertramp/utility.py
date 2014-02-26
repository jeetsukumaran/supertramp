#! /usr/bin/env python

##############################################################################
##
##  Copyright 2010-2014 Jeet Sukumaran.
##  All rights reserved.
##
##  Redistribution and use in source and binary forms, with or without
##  modification, are permitted provided that the following conditions are met:
##
##      * Redistributions of source code must retain the above copyright
##        notice, this list of conditions and the following disclaimer.
##      * Redistributions in binary form must reproduce the above copyright
##        notice, this list of conditions and the following disclaimer in the
##        documentation and/or other materials provided with the distribution.
##      * The names of its contributors may not be used to endorse or promote
##        products derived from this software without specific prior written
##        permission.
##
##  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
##  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
##  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
##  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN OR MARK T. HOLDER
##  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
##  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
##  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
##  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
##  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
##  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
##  POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

from io import StringIO
import logging

_LOGGING_LEVEL_ENVAR = "SUPERTRAMP_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR = "SUPERTRAMP_LOGGING_FORMAT"

class RunLogger(object):

    def __init__(self, **kwargs):
        self.name = kwargs.get("name", "RunLog")
        self._log = logging.getLogger(self.name)
        self._log.setLevel(logging.DEBUG)
        if kwargs.get("log_to_stderr", True):
            ch1 = logging.StreamHandler()
            stderr_logging_level = self.get_logging_level(kwargs.get("stderr_logging_level", logging.INFO))
            ch1.setLevel(stderr_logging_level)
            ch1.setFormatter(self.get_default_formatter())
            self._log.addHandler(ch1)
        if kwargs.get("log_to_file", True):
            log_stream = kwargs.get("log_stream", \
                open(kwargs.get("log_path", self.name + ".log"), "w"))
            ch2 = logging.StreamHandler(log_stream)
            file_logging_level = self.get_logging_level(kwargs.get("file_logging_level", logging.DEBUG))
            ch2.setLevel(file_logging_level)
            ch2.setFormatter(self.get_default_formatter())
            self._log.addHandler(ch2)

    def get_logging_level(self, level=None):
        if level in [logging.NOTSET, logging.DEBUG, logging.INFO, logging.WARNING,
            logging.ERROR, logging.CRITICAL]:
            return level
        elif level is not None:
            level_name = str(level).upper()
        elif _LOGGING_LEVEL_ENVAR in os.environ:
            level_name = os.environ[_LOGGING_LEVEL_ENVAR].upper()
        else:
            level_name = "NOTSET"
        if level_name == "NOTSET":
            level = logging.NOTSET
        elif level_name == "DEBUG":
            level = logging.DEBUG
        elif level_name == "INFO":
            level = logging.INFO
        elif level_name == "WARNING":
            level = logging.WARNING
        elif level_name == "ERROR":
            level = logging.ERROR
        elif level_name == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
        return level

    def get_default_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_rich_formatter(self):
        f = logging.Formatter("[%(asctime)s] %(filename)s (%(lineno)d): %(levelname) 8s: %(message)s")
        f.datefmt='%Y-%m-%d %H:%M:%S'
        return f

    def get_simple_formatter(self):
        return logging.Formatter("%(levelname) 8s: %(message)s")

    def get_raw_formatter(self):
        return logging.Formatter("%(message)s")

    def get_logging_formatter(self, format=None):
        if format is not None:
            format = format.upper()
        elif _LOGGING_FORMAT_ENVAR in os.environ:
            format = os.environ[_LOGGING_FORMAT_ENVAR].upper()
        if format == "RICH":
            logging_formatter = self.get_rich_formatter()
        elif format == "SIMPLE":
            logging_formatter = self.get_simple_formatter()
        elif format == "NONE":
            logging_formatter = self.get_raw_formatter()
        else:
            logging_formatter = self.get_default_formatter()
        if logging_formatter is not None:
            logging_formatter.datefmt='%H:%M:%S'

    def debug(self, msg, *args, **kwargs):
        self._log.debug(msg, *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        self._log.info(msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        self._log.warning(msg, *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self._log.error(msg, *args, **kwargs)

    def critical(self, msg, *args, **kwargs):
        self._log.critical(msg, *args, **kwargs)

