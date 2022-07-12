# -*- coding: utf-8 -*-

from __future__ import absolute_import

# set the warning format to be on a single line
import sys
import logging
import warnings as _warnings
from os import name as _name
from os.path import abspath as _abspath
from os.path import dirname as _dirname

logging_hash = {
    "debug":logging.DEBUG,
    "critical":logging.CRITICAL,
    "error":logging.ERROR,
    "warning":logging.WARNING,
    "info":logging.INFO
}

#Configuing modelseedpy logger
#logger = logging.getLogger(__name__)

from ModelSEEDReconstruction.basemodule import SDKSupportModule
