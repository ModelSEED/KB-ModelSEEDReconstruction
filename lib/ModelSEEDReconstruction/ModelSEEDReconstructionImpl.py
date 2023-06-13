# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import sys
sys.path.append("/deps/KBBaseModules/")
import json
from os.path import exists
from ModelSEEDReconstruction.modelseedrecon import ModelSEEDRecon
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil

logger = logging.getLogger(__name__)
logger.setLevel(
    logging.INFO
) 
#END_HEADER


class ModelSEEDReconstruction:
    '''
    Module Name:
    ModelSEEDReconstruction

    Module Description:
    A KBase module: ModelSEEDReconstruction
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    def build_report(self,output):
        report_shock_id = self.dfu.file_to_shock({'file_path': output["file_path"],'pack': 'zip'})['shock_id']
        output["report_params"]["html_links"][0]["shock_id"] = report_shock_id
        repout = self.kbreport.create_extended_report(output["report_params"])
        return {"report_name":output["report_params"]["report_object_name"],"report_ref":repout["ref"],'workspace_name':self.msrecon.ws_name}
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config["ATP_media_workspace"] = "94026"
        if "appdev" in self.config["workspace-url"]:
            self.config["ATP_media_workspace"] = "68393" 
        config["version"] = self.VERSION
        self.msrecon = ModelSEEDRecon(config,"/kb/module",config['scratch'],os.environ['KB_AUTH_TOKEN'],{},os.environ['SDK_CALLBACK_URL'])
        #END_CONSTRUCTOR
        pass

    def build_metabolic_models(self, ctx, params):
        """
        Function builds metabolic models from input genomes
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String,
           parameter "workspace" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN build_metabolic_models
        #Processing parameters
        output = self.msrecon.build_metabolic_models(params)
        output = self.build_report(output)
        #END build_metabolic_models

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method build_metabolic_models return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def gapfill_metabolic_models(self, ctx, params):
        """
        Function gapfills metabolic models from input genomes
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String,
           parameter "workspace" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN gapfill_metabolic_models
        output = self.msrecon.gapfill_metabolic_models(params)
        output = self.build_report(output)
        #END gapfill_metabolic_models

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method gapfill_metabolic_models return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
