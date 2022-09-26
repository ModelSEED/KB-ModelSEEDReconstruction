# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

import sys
import json
import cobra
import re
import uuid
from os.path import exists
import cobrakbase
from optlang.symbolics import Zero, add
from modelseedpy import MSPackageManager,MSGenome, MSMedia, MSModelUtil, MSBuilder, MSGapfill, FBAHelper, MSGrowthPhenotypes, MSModelUtil, MSATPCorrection
from cobrakbase.core.kbasefba.newmodeltemplate_builder import NewModelTemplateBuilder
from cobrakbase.core.kbasefba.fbamodel_from_cobra import CobraModelConverter
from modelseedpy.helpers import get_template
from modelseedpy.core.msgenomeclassifier import MSGenomeClassifier
from modelseedpy.core.mstemplate import MSTemplateBuilder
from ModelSEEDReconstruction.sdkhelper import SDKHelper
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
import pickle
import pandas as pd

def fix_genomescale_template(gs_template,core_template):
    for cpd in core_template.compcompounds:
        if cpd.id not in gs_template.compcompounds:
            gs_template.compcompounds.append(cpd)
    for rxn in core_template.reactions:
        if rxn.id in gs_template.reactions:
            gs_template.reactions._replace_on_id(rxn)
        else:
            gs_template.reactions.append(rxn)
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
    def save_model(self,mdlutl,workspace):
        fbamodel = CobraModelConverter(mdlutl.model, mdlutl.model.genome, mdlutl.model.template).build()  # later assign core template to None
        json = fbamodel.get_data()
        mdlutl.create_kb_gapfilling_data(json)
        params = {
            'objects': [{
                'data': json,
                'name': mdlutl.model.id,
                'type': "KBaseFBA.FBAModel",
                'meta': {},
                'provenance': [{
                    'description': '',
                    'input_ws_objects': self.input_objects,
                    'method': self.method,
                    'script_command_line': "",
                    'method_params': [self.params],
                    'service': 'ModelSEEDReconstruction',
                    'service_ver': ModelSEEDReconstruction.VERSION,
                    # 'time': '2015-12-15T22:58:55+0000'
                }]
            }]
        }
        if isinstance(workspace, int):
            params["id"] = workspace
        else:
            params["workspace"] = workspace
        self.kbase_api.ws_client.save_objects(params)
        self.obj_created.append(SDKHelper.create_ref(mdlutl.model.id,workspace))
    
    def build_report(self,table,workspace):
        #columns=column_list
        html_data = f"""
    <html>
    <header>
        <link href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.min.css" rel="stylesheet">
    </header>
    <body>
    {table.to_html(escape=False,notebook=False,table_id="table",index=False,justify="left")}
    <script src="https://code.jquery.com/jquery-3.6.0.slim.min.js" integrity="sha256-u7e5khyithlIdTpu22PHhENmPcRdFiHRjhAuHcs05RI=" crossorigin="anonymous"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
    <script>
        $(document).ready( function () {{
            $('#table').DataTable({{
                // paging: false,    
                // scrollY: 400,
            }});
        }});
    </script>
    </body>
    </html>
    """
        report_name = str(uuid.uuid4())
        if not isinstance(workspace, str):
            workspace = str(workspace)
        html_report_folder = os.path.join(self.shared_folder, 'htmlreport')
        os.makedirs(html_report_folder, exist_ok=True)
        with open(os.path.join(html_report_folder, 'index.html'), 'w') as f:
            f.write(html_data)
        report_shock_id = self.dfu.file_to_shock({'file_path': html_report_folder,'pack': 'zip'})['shock_id']
        html_output = {
            'name' : 'index.html',
            'shock_id': report_shock_id
        }
        report_params = {
            'objects_created': [self.obj_created],
            'workspace_name': workspace,
            'html_links': [{
                'name' : 'index.html',
                'shock_id': report_shock_id
            }],
            'direct_html_link_index': 0,
            'html_window_height': 700,
            'report_object_name': report_name
        }
        report = KBaseReport(self.callback_url, token=self.token)
        repout = report.create_extended_report(report_params)
        return {"report_name":report_name,"report_ref":repout["ref"],'workspace_name':workspace}
    
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.token = os.environ['KB_AUTH_TOKEN']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.dfu = DataFileUtil(self.callback_url)
        self.kbase_api = cobrakbase.KBaseAPI(token=self.token,config=self.config)
        self.obj_created = []
        self.input_objects = []
        self.method = None
        self.params = None
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
        if self.method == None:
            self.method = "build_metabolic_models"
            self.params = params
        SDKHelper.validate_args(params,["workspace"],{
            "genome_refs":[],
            "run_gapfilling":False,
            "atp_safe":True,
            "forced_atp_list":[],
            "gapfilling_media_list":[],
            "suffix":".mdl",
            "core_template":"auto",
            "gs_template":"auto",
            "gs_template_ref":None,
            "template_reactions_only":True,
            "output_core_models":False
        })
        #Preloading core and preselected template
        templates = {
            "core" : self.kbase_api.get_from_ws("core_template_sulfur3","NewKBaseModelTemplates"),
            "gp" : None,
            "gn" : None,
            "custom": None
        }
        if params["gs_template"] == "custom":
            templates["custom"] = self.kbase_api.get_from_ws(params["gs_template_ref"],None)
        #Initializing classifier
        cls_pickle = "/kb/module/data/knn_ACNP_RAST_filter.pickle"
        cls_features = "/kb/module/data/knn_ACNP_RAST_filter_features.json"
        with open(cls_pickle, 'rb') as fh:
            model_filter = pickle.load(fh)
        with open(cls_features, 'r') as fh:
            features = json.load(fh)
        genome_classifier = MSGenomeClassifier(model_filter, features)
        #Initializing output data tables
        result_table = pd.DataFrame({})
        default_output = {"Model":None,"Genome":None,"Genes":None,"Class":None,
                          "Model genes":None,"Reactions":None,"ATP yeilds":None,
                          "Core GF":None,"GS GF":None,"Auxotrophy":None,"Growth":None,"Comments":None}
        #Retrieving genomes and building models one by one
        for i,gen_ref in enumerate(params["genome_refs"]):
            self.input_objects.append(gen_ref)
            genome = self.kbase_api.get_from_ws(gen_ref)
            #Initializing output row
            current_output = default_output.copy()
            comments = []
            gid = genome.id
            current_output["Model"] = gid+params["suffix"]
            current_output["Genome"] = genome.info.metadata["Name"]
            current_output["Genes"] = genome.info.metadata["Number of Protein Encoding Genes"]
            #Pulling annotation priority
            comments.append("Other annotation priorities not supported by this app yet. Using RAST.")
            #Retrieving genome annotations
            template_type = params["gs_template"]
            if template_type == "auto":
                current_output["Class"] = genome_classifier.classify(genome)
                if current_output["Class"] == "P":
                    template_type = "gp"
                elif current_output["Class"] == "N" or current_output["Class"] == "--":
                    template_type = "gn"
                if template_type not in templates:
                    current_output["Comments"] = "Template type "+template_type+" not recognized"
                    result_table.append(current_output)
                    next
                elif templates[template_type] == None:
                    if template_type == "gn":
                        templates[template_type] = self.kbase_api.get_from_ws("GramNegModelTemplateV4","NewKBaseModelTemplates")
                    if template_type == "gp":
                        templates[template_type] = self.kbase_api.get_from_ws("GramPosModelTemplateV4","NewKBaseModelTemplates")
                    fix_genomescale_template(templates[template_type],templates["core"])
            curr_template = templates[template_type]
            #Building model
            mdllist = []
            mdllist.append(MSBuilder(genome, curr_template).build(gid+params["suffix"], '0', False, False))
            mdllist[0].template = curr_template
            mdllist.append(MSBuilder(genome, templates["core"]).build(gid+".core"+params["suffix"], '0', False, False))
            mdllist[1].template = templates["core"]
            for i,mdl in enumerate(mdllist):
                mdl.genome = genome
                mdlutl = MSModelUtil(mdl)
                mdllist[i] = mdlutl
                if params["atp_safe"] and mdlutl == mdllist[0]:
                    mdlutl.set_atputl(MSATPCorrection.build_default(mdlutl,templates["core"],atp_medias=[],forced_media=params["forced_atp_list"],default_media_path="/kb/module/data/atp_medias.tsv"))
                    mdlutl.atputl.evaluate_growth_media()
                    mdlutl.atputl.determine_growth_media()
                    mdlutl.atputl.apply_growth_media_gapfilling()
                    mdlutl.atputl.evaluate_growth_media()
                    mdlutl.atputl.expand_model_to_genome_scale()
                    if mdlutl == mdllist[0]:
                        tests = mdlutl.atputl.build_tests()
                        current_output["ATP yeilds"] = ""
                        for test in tests:
                            if len(current_output["ATP yeilds"]) > 0:
                                current_output["ATP yeilds"] += "; "
                            current_output["ATP yeilds"] += test["media"].id+":"+str(test["threshold"])
                    else:
                        current_output["Core GF"] = mdlutl.gapfilled_reaction_count()
                elif mdlutl == mdllist[0]:
                    current_output["ATP yeilds"] = "NA"
                    current_output["Core GF"] = "NA"    
            #Running gapfilling
            current_output["GS GF"] = "NA"
            current_output["Auxotrophy"] = "NA"
            if params["run_gapfilling"]:
                self.gapfill_metabolic_models(ctx,{
                    "media_list":params["gapfilling_media_list"],#
                    "model_objs":[mdllist[0]],#
                    "atp_safe":params["atp_safe"],#
                    "workspace":params["workspace"],#
                    "suffix":params["suffix"],#
                    "default_objective":"bio1",#
                    "output_data":{mdllist[0]:current_output},#
                    "forced_atp_list":params["forced_atp_list"],
                    "templates":[curr_template] ,
                    "internal_call":True
                })
            else:
                self.save_model(mdllist[0],params["workspace"])
            if params["output_core_models"]:
                self.save_model(mdllist[1],params["workspace"])
            #Filling in model output
            current_output["Reactions"] = len(mdllist[0].model.reactions)
            current_output["Model genes"] = len(mdllist[0].model.genes)
            mdllist[0].model.objective = "bio1"
            current_output["Growth"] = mdllist[0].model.slim_optimize()
            result_table = result_table.append(current_output, ignore_index = True)
        output = self.build_report(result_table,params["workspace"])
        output["data"] = result_table.to_json()
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
        #Processing parameters
        if self.method == None:
            self.method = "gapfill_metabolic_models"
            self.params = params
        SDKHelper.validate_args(params,["media_list","workspace"],{
            "model_list":None,
            "model_objectives":[],
            "model_objs":[],
            "atp_safe":True,
            "suffix":".gf",
            "forced_atp_list":[],
            "templates":None,
            "source_models":[],
            "limit_medias":[],
            "limit_objectives":[],
            "limit_thresholds":[],
            "is_max_limits":[],
            "minimum_objective":0.01,
            "reaction_exlusion_list":[],
            "default_objective":"bio1",
            "kbmodel_hash":{},
            "output_data":None,
            "internal_call":False
        })
        result_table = pd.DataFrame({})
        default_output = {"Model":None,"Genome":None,"Genes":None,"Class":None,
            "Model genes":None,"Reactions":None,"ATP yeilds":None,
            "Core GF":None,"GS GF":None,"Auxotrophy":None,"Growth":None,"Comments":None}
        #Retrieving models if not provided already
        if "model_objs" not in params or len(params["model_objs"]) == 0:
            params["model_objs"] = []
            for mdl_ref in params["model_list"]:
                self.input_objects.append(mdl_ref)
                kbmodel = self.kbase_api.get_object(mdl_ref,None)
                model = self.kbase_api.get_from_ws(mdl_ref,None)
                model.genome = self.kbase_api.get_from_ws(kbmodel["genome_ref"],None)
                model.template = self.kbase_api.get_from_ws(kbmodel["template_ref"],None)
                model.id = model.info[0] + params["suffix"]
                params["model_objs"].append(MSModelUtil(model))
        #Retrieving media objects from references
        media_objs = []
        for media_ref in params["media_list"]:
            self.input_objects.append(media_ref)
            media = self.kbase_api.get_from_ws(media_ref,None)
            media_objs.append(media)
        #Compiling additional tests
        additional_tests = []
        for i,limit_media in enumerate(params["limit_medias"]):
            additional_tests.append({
                "objective":params["limit_objectives"][i],
                "media":limit_media,
                "is_max_threshold":params["is_max_limits"][i],
                "threshold":params["limit_thresholds"][i]
            })
        #Iterating over each model and running gapfilling
        for i,mdlutl in enumerate(params["model_objs"]):
            #Setting the objective
            if i < len(params["model_objectives"]):
                if not params["model_objectives"][i]:
                    params["model_objectives"][i] = params["default_objective"]
            else:
                params["model_objectives"].append(params["default_objective"])
            #Creating gapfilling object and configuring solver
            #mdlutl.model.solver = config["solver"]
            mdlutl.set_gfutl(MSGapfill.build_default(mdlutl,params["atp_safe"],params["templates"],params["source_models"],
                     additional_tests,params["reaction_exlusion_list"]))
            #Iterating over all media specified for gapfilling
            for media in media_objs:
                #Gapfilling
                gfresults = mdlutl.gfutl.run_gapfilling(media,params["model_objectives"][i],
                    params["minimum_objective"])
                mdlutl.gfutl.integrate_gapfill_solution(gfresults)
                mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
                #solution = mdlutl.model.optimize()
            #Saving completely gapfilled model
            self.save_model(mdlutl,params["workspace"])
        if not params["internal_call"]:
            output = self.build_report(result_table,params["workspace"])
            output["data"] = result_table.to_json()
        else:
            output["data"] = result_table
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
