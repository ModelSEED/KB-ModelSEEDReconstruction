from __future__ import absolute_import

import os
import sys
import uuid
import logging
import json
import pandas as pd
import cobrakbase
from optlang.symbolics import Zero, add
from modelseedpy import MSPackageManager,MSGenome, MSMedia, MSModelUtil, MSBuilder, MSGapfill, FBAHelper, MSGrowthPhenotypes, MSModelUtil, MSATPCorrection
from cobrakbase.core.kbasefba.newmodeltemplate_builder import NewModelTemplateBuilder
from modelseedpy.helpers import get_template
from modelseedpy.core.msgenomeclassifier import MSGenomeClassifier
from modelseedpy.core.mstemplate import MSTemplateBuilder
from ModelSEEDReconstruction.basemodelingmodule import BaseModelingModule
import pickle

logger = logging.getLogger(__name__)

def fix_genomescale_template(gs_template,core_template):
    for cpd in core_template.compcompounds:
        if cpd.id not in gs_template.compcompounds:
            gs_template.compcompounds.append(cpd)
    for rxn in core_template.reactions:
        if rxn.id in gs_template.reactions:
            gs_template.reactions._replace_on_id(rxn)
        else:
            gs_template.reactions.append(rxn)

class ModelSEEDRecon(BaseModelingModule):
    def __init__(self,name,ws_client,working_dir,module_dir,config):
        BaseModelingModule.__init__(self,name,ws_client,working_dir,config)
        self.module_dir = module_dir
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        
    def build_metabolic_models(self,params):
        self.initialize_call("build_metabolic_models",params,True)
        if "genomes" in params:
            params["genome_refs"] = params["genomes"]["genome_refs"]
        if "reconstruction_parameters" in params:
            params["suffix"] = params["reconstruction_parameters"]["suffix"]
            params["output_core_models"] = params["reconstruction_parameters"]["output_core_models"]
        if "gapfilling_parameters" in params:
            params["run_gapfilling"] = params["gapfilling_parameters"]["run_gapfilling"]
            params["gapfilling_media_list"] = params["gapfilling_parameters"]["gapfilling_media_list"]
        if "atp_safe_parameters" in params:
            params["atp_safe"] = params["atp_safe_parameters"]["atp_safe"]
            params["forced_atp_list"] = params["atp_safe_parameters"]["forced_atp_list"]
            params["automated_atp_evaluation"] = params["atp_safe_parameters"]["automated_atp_evaluation"]
        if "template_parameters" in params:
            params["gs_template"] = params["template_parameters"]["gs_template"]
            params["gs_template_ref"] = params["template_parameters"]["gs_template_ref"]
            params["core_template_ref"] = params["template_parameters"]["core_template_ref"]
        self.validate_args(params,["workspace"],{
            "genome_refs":[],
            "run_gapfilling":False,
            "atp_safe":True,
            "forced_atp_list":[],
            "gapfilling_media_list":[],
            "suffix":".mdl",
            "core_template":"auto",
            "gs_template":"auto",
            "gs_template_ref":None,
            "core_template_ref":None,
            "template_reactions_only":True,
            "output_core_models":False,
            "automated_atp_evaluation":True,
            "atp_medias":[],
            "load_default_medias":True,
            "max_gapfilling":10,
            "gapfilling_delta":0
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
        cls_pickle = self.module_dir+"/data/knn_ACNP_RAST_filter.pickle"
        cls_features = self.module_dir+"/data/knn_ACNP_RAST_filter_features.json"
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
                    fix_genomescale_template(templates[template_type],templates["core"])#Move to MSTemplate?
            curr_template = templates[template_type]
            #Building model
            mdl = MSBuilder(genome, curr_template).build(gid+params["suffix"], '0', False, False)
            mdl.genome = genome#Move into MSBuilder?
            mdl.template = curr_template#Move into MSBuilder?
            current_output["ATP yeilds"] = "NA"
            current_output["Core GF"] = "NA"  
            mdlutl = MSModelUtil.get(mdl)
            if params["atp_safe"]:
                atpcorrection = MSATPCorrection(mdlutl,templates["core"],params["atp_medias"],load_default_medias=params["load_default_medias"],max_gapfilling=params["max_gapfilling"],gapfilling_delta=params["gapfilling_delta"],forced_media=params["forced_atp_list"],default_media_path=self.module_dir+"/data/atp_medias.tsv")
                tests = atpcorrection.run_atp_correction()
                current_output["ATP yeilds"] = ""
                for test in tests:
                    if len(current_output["ATP yeilds"]) > 0:
                        current_output["ATP yeilds"] += "; "
                    current_output["ATP yeilds"] += test["media"].id+":"+str(test["threshold"])
                current_output["Core GF"] = len(atpcorrection.cumulative_core_gapfilling)
            #Running gapfilling
            current_output["GS GF"] = "NA"
            current_output["Auxotrophy"] = "NA"
            if params["run_gapfilling"]:
                self.gapfill_metabolic_models({
                    "media_list":params["gapfilling_media_list"],#
                    "model_objs":[mdlutl],#
                    "atp_safe":params["atp_safe"],#
                    "workspace":params["workspace"],#
                    "suffix":params["suffix"],#
                    "default_objective":"bio1",#
                    "output_data":{mdlutl:current_output},#
                    "forced_atp_list":params["forced_atp_list"],
                    "templates":[curr_template] ,
                    "internal_call":True
                })
                current_output["GS GF"] = len(mdlutl.gfutl.cumulative_gapfilling)
            else:
                self.save_model(mdlutl,params["workspace"])
            current_output["Reactions"] = len(mdlutl.model.reactions)
            current_output["Model genes"] = len(mdlutl.model.genes)
            mdlutl.model.objective = "bio1"
            current_output["Growth"] = mdlutl.model.slim_optimize()
            #if params["output_core_models"]:
            #TODO: Remove noncore reactions and change biomass and change model ID and then resave
            #Filling in model output
            result_table = result_table.append(current_output, ignore_index = True)
        output = self.build_report(result_table)
        output["data"] = result_table.to_json()
    
    def gapfill_metabolic_models(self,params):
        self.initialize_call("gapfill_metabolic_models",params,True)
        self.validate_args(params,["media_list","workspace"],{
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
            msgapfill = MSGapfill(mdlutl,params["templates"],params["source_models"],
                     additional_tests,blacklist=params["reaction_exlusion_list"])
            #Iterating over all media specified for gapfilling
            mdlutl.gfutl.cumulative_gapfilling = []
            for media in media_objs:
                #Gapfilling
                gfresults = msgapfill.run_gapfilling(media,params["model_objectives"][i],
                    params["minimum_objective"])
                msgapfill.integrate_gapfill_solution(gfresults,mdlutl.gfutl.cumulative_gapfilling)
                #mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
                #solution = mdlutl.model.optimize()
            #Saving completely gapfilled model
            self.save_model(mdlutl,params["workspace"])
        output = {}
        if not params["internal_call"]:
            output = self.build_report(result_table)
            output["data"] = result_table.to_json()
            
    def build_report(self,table):
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
        html_report_folder = os.path.join(self.working_dir, 'htmlreport')
        os.makedirs(html_report_folder, exist_ok=True)
        with open(os.path.join(html_report_folder, 'index.html'), 'w') as f:
            f.write(html_data)
        return {
            'data':table,
            'file_path':os.path.join(html_report_folder, 'index.html'),
            'report_params':{
                'objects_created': self.obj_created,
                'workspace_name': self.ws_name,
                'html_links': [{
                    'name' : 'index.html',
                    'shock_id': None
                }],
                'direct_html_link_index': 0,
                'html_window_height': 700,
                'report_object_name': report_name
            }
        }