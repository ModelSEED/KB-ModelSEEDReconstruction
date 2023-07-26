from __future__ import absolute_import

import os
import sys
import uuid
import logging
import json
import jinja2
import pandas as pd
from optlang.symbolics import Zero, add
from cobrakbase.core.kbasefba import FBAModel
from modelseedpy import MSPackageManager,MSGenome, MSMedia, MSModelUtil, MSBuilder, MSGapfill, FBAHelper, MSGrowthPhenotypes, MSModelUtil, MSATPCorrection,MSModelReport
from modelseedpy.helpers import get_template
from modelseedpy.core.msgenomeclassifier import MSGenomeClassifier
from modelseedpy.core.mstemplate import MSTemplateBuilder
from kbbasemodules.basemodelingmodule import BaseModelingModule

logger = logging.getLogger(__name__)

class ModelSEEDRecon(BaseModelingModule):
    def __init__(self,config,module_dir="/kb/module",working_dir=None,token=None,clients={},callback=None):
        BaseModelingModule.__init__(self,"ModelSEEDReconstruction",config,module_dir,working_dir,token,clients,callback)
        self.module_dir = module_dir
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        
    def build_metabolic_models(self,params):
        self.initialize_call("build_metabolic_models",params,True)
        if "genomes" in params:
            params["genome_refs"] = params["genomes"]["genome_refs"]
        if "reconstruction_parameters" in params:
            params["suffix"] = params["reconstruction_parameters"]["suffix"]
            #params["output_core_models"] = params["reconstruction_parameters"]["output_core_models"]
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
            "gapfilling_media_list":["KBaseMedia/Complete"],
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
            "gapfilling_delta":0,
            "return_model_objects":False,
            "return_data":False,
            "save_report_to_kbase":True
        })
        #Preloading core and preselected template
        template_type = params["gs_template"]
        if params["gs_template_ref"]:
            self.templates["custom"] = self.get_template(params["gs_template_ref"],None)
            template_type = "custom"
        if params["core_template_ref"]:
            self.templates["core"] = self.get_template(params["core_template_ref"],None)
        #Initializing classifier
        genome_classifier = self.get_classifier()
        #Initializing output data tables
        result_table = pd.DataFrame({})
        default_output = {"Model":None,"Genome":None,"Genes":None,"Class":None,
                          "Model genes":None,"Reactions":None,"ATP yeilds":None,
                          "Core GF":None,"GS GF":None,"Auxotrophy":None,"Growth":None,"Comments":None}
        #Retrieving genomes and building models one by one
        mdllist = []
        for i,gen_ref in enumerate(params["genome_refs"]):
            genome = self.get_msgenome(gen_ref)
            #Initializing output row
            current_output = default_output.copy()
            comments = []
            gid = genome.id
            current_output["Model"] = gid+params["suffix"]
            current_output["Genome"] = genome.info.metadata["Name"]
            current_output["Genes"] = genome.info.metadata["Number of Protein Encoding Genes"]
            #Pulling annotation priority
            comments.append("Other annotation priorities not supported by this app yet. Using RAST.")
            if template_type == "auto":
                current_output["Class"] = genome_classifier.classify(genome)
                print(current_output["Class"])
                if current_output["Class"] == "P":
                    template_type = "gp"
                elif current_output["Class"] == "N" or current_output["Class"] == "--":
                    template_type = "gn"
                elif current_output["Class"] == "A":
                    template_type = "ar"
                if template_type not in self.templates:
                    current_output["Comments"] = "Template type "+template_type+" not recognized"
                    result_table = result_table.append(current_output, ignore_index = True)
                    next
                elif self.templates[template_type] == None:
                    if template_type == "gn":
                        self.templates[template_type] = self.get_template("GramNegModelTemplateV5","NewKBaseModelTemplates")#,templates["core"])
                    if template_type == "gp":
                        self.templates[template_type] = self.get_template("GramPosModelTemplateV5","NewKBaseModelTemplates")#,templates["core"])
            curr_template = self.templates[template_type]
            #Building model            
            base_model = FBAModel({'id':current_output["Model"], 'name':genome.scientific_name})
            mdl = MSBuilder(genome, curr_template).build(base_model, '0', False, False)
            mdl.genome = genome
            mdl.template = curr_template
            mdl.core_template_ref = str(self.templates["core"].info)
            mdl.genome_ref = str(genome.info)
            mdl.template_ref = str(curr_template.info)
            current_output["ATP yeilds"] = "NA"
            current_output["Core GF"] = "NA" 
            mdlutl = MSModelUtil.get(mdl)
            if params["atp_safe"]:
                atpcorrection = MSATPCorrection(mdlutl,self.templates["core"],params["atp_medias"],load_default_medias=params["load_default_medias"],max_gapfilling=params["max_gapfilling"],gapfilling_delta=params["gapfilling_delta"],forced_media=params["forced_atp_list"],default_media_path=self.module_dir+"/data/atp_medias.tsv")
                tests = atpcorrection.run_atp_correction()
                current_output["ATP yeilds"] = ""
                for test in tests:
                    if len(current_output["ATP yeilds"]) > 0:
                        current_output["ATP yeilds"] += "; "
                    current_output["ATP yeilds"] += test["media"].id+":"+str(test["threshold"])
                current_output["Core GF"] = len(atpcorrection.cumulative_core_gapfilling)
            #Setting the model ID so the model is saved with the correct name in KBase
            mdlutl.wsid = gid+params["suffix"]
            #Running gapfilling
            current_output["GS GF"] = "NA"
            current_output["Auxotrophy"] = "NA"
            if params["run_gapfilling"]:
                self.gapfill_metabolic_models({
                    "media_list":params["gapfilling_media_list"],#
                    "model_objs":[mdlutl],#
                    "atp_safe":params["atp_safe"],#
                    "workspace":params["workspace"],#
                    "suffix":"",#
                    "default_objective":"bio1",#
                    "output_data":{mdlutl:current_output},#
                    "forced_atp_list":params["forced_atp_list"],
                    "templates":[curr_template],
                    "internal_call":True
                })
            else:
                self.save_model(mdlutl,params["workspace"],None)
                mdlutl.model.objective = "bio1"
                mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(None)
                current_output["Growth"] = "Complete:"+str(mdlutl.model.slim_optimize())
            current_output["Reactions"] = mdlutl.nonexchange_reaction_count()
            current_output["Model genes"] = len(mdlutl.model.genes)
            #if params["output_core_models"]:
            #TODO: Remove noncore reactions and change biomass and change model ID and then resave
            #Filling in model output
            result_table = result_table.append(current_output, ignore_index = True)
            mdllist.append(mdlutl)
        output = {}
        if params["save_report_to_kbase"]:
            self.build_dataframe_report(result_table,mdllist)
            output = self.save_report_to_kbase()
        if params["return_data"]:
            output["data"] = result_table.to_json()
        if params["return_model_objects"]:
            output["model_objs"] = mdllist
        return output
    
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
            "internal_call":False,
            "atp_medias":[],
            "load_default_atp_medias":True,
            "max_atp_gapfilling":0,
            "gapfilling_delta":0,
            "return_data":False,
            "save_report_to_kbase":True
        })
        result_table = pd.DataFrame({})
        default_output = {"Model":None,"Genome":None,"Genes":None,"Class":None,
            "Model genes":None,"Reactions":None,"ATP yeilds":None,
            "Core GF":None,"GS GF":None,"Auxotrophy":None,"Growth":None,"Comments":None}
        #Retrieving models if not provided already
        if "model_objs" not in params or len(params["model_objs"]) == 0:
            params["model_objs"] = []
            for mdl_ref in params["model_list"]:
                params["model_objs"].append(self.get_model(mdl_ref))
        #Retrieving media objects from references
        params["media_objs"] = []
        for media_ref in params["media_list"]:
            self.input_objects.append(media_ref)
            media = self.get_media(media_ref,None)
            params["media_objs"].append(media)
        #Compiling additional tests
        additional_tests = []
        for i,limit_media in enumerate(params["limit_medias"]):
            additional_tests.append({
                "objective":params["limit_objectives"][i],
                "media":self.get_media(limit_media,None),
                "is_max_threshold":params["is_max_limits"][i],
                "threshold":params["limit_thresholds"][i]
            })
        #Iterating over each model and running gapfilling
        for i,mdlutl in enumerate(params["model_objs"]):
            current_output = default_output.copy()
            current_output["Model"] = mdlutl.wsid+params["suffix"]    
            if params["output_data"] and mdlutl in params["output_data"]:
                current_output = params["output_data"][mdlutl]
            #Setting the objective
            if i < len(params["model_objectives"]):
                if not params["model_objectives"][i]:
                    params["model_objectives"][i] = params["default_objective"]
            else:
                params["model_objectives"].append(params["default_objective"])
            #Computing tests for ATP safe gapfilling
            if params["atp_safe"]:
                if not mdlutl.atputl:
                    atpcorrection = MSATPCorrection(mdlutl,self.templates["core"],params["atp_medias"],load_default_medias=params["load_default_atp_medias"],max_gapfilling=params["max_atp_gapfilling"],gapfilling_delta=params["gapfilling_delta"],forced_media=params["forced_atp_list"],default_media_path=self.module_dir+"/data/atp_medias.tsv")
                    tests = atpcorrection.run_atp_correction()
                    additional_tests.extend(tests)
                else:
                    tests = mdlutl.atputl.build_tests()
                    additional_tests.extend(tests)
            #Creating gapfilling object and configuring solver
            #mdlutl.model.solver = config["solver"]
            if not params["templates"]:
                params["templates"] = [self.get_template(mdlutl.model.template_ref)]
            msgapfill = MSGapfill(mdlutl,params["templates"],params["source_models"],
                     additional_tests,blacklist=params["reaction_exlusion_list"],default_target=params["model_objectives"][i],minimum_obj=params["minimum_objective"])
            #Running gapfilling in all conditions
            mdlutl.gfutl.cumulative_gapfilling = []
            growth_array = []
            solutions = msgapfill.run_multi_gapfill(params["media_objs"],default_minimum_objective=params["minimum_objective"])
            for media in params["media_objs"]:
                if media in solutions:
                    msgapfill.integrate_gapfill_solution(solutions[media],mdlutl.gfutl.cumulative_gapfilling)
                mdlutl.model.objective = params["model_objectives"][i]
                mdlutl.pkgmgr.getpkg("KBaseMediaPkg").build_package(media)
                growth_array.append(media.id+":"+str(mdlutl.model.slim_optimize()))
            current_output["Growth"] = "<br>".join(growth_array)
            current_output["GS GF"] = len(mdlutl.gfutl.cumulative_gapfilling)
            current_output["Reactions"] = mdlutl.nonexchange_reaction_count()
            current_output["Model genes"] = len(mdlutl.model.genes)
            #current_output["GF reasons"] = ""
            #mdlutl.gfutl.link_gapfilling_to_biomass()
            #for item in mdlutl.gfutl.cumulative_gapfilling:
            #    if len(item) > 2 and len(item[2]) > 0:
            #        if len(current_output["GF reasons"]) > 0:
            #            current_output["GF reasons"] += "<br>"
            #        current_output["GF reasons"] += item[0].id+item[1]+":"+",".join(item[2])
            #Saving completely gapfilled model
            self.save_model(mdlutl,params["workspace"],None,params["suffix"])
            if not params["internal_call"]:
                result_table = result_table.append(current_output, ignore_index = True)
        output = {}
        if not params["internal_call"]:
            if params["save_report_to_kbase"]:
                self.build_dataframe_report(result_table,params["model_objs"])
                output = self.save_report_to_kbase()
            if params["return_data"]:
                output["data"] = result_table.to_json()
        return output
            
    def build_dataframe_report(self,table,model_objs=None):        
        context = {
            "initial_model":table.iloc[0]["Model"]
        }
        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(self.module_dir+"/data/"),
            autoescape=jinja2.select_autoescape(['html', 'xml']))
        html = env.get_template("ReportTemplate.html").render(context)
        os.makedirs(self.working_dir+"/html", exist_ok=True)
        if model_objs:
            msmodrep = MSModelReport()
            for model in model_objs:
                msmodrep.build_report(model,self.working_dir+"/html/"+model.wsid+".html")  
        with open(self.working_dir+"/html/index.html", 'w') as f:
            f.write(html)
        #Creating data table file
        for index, row in table.iterrows():
            table.at[index,'Model'] = '<a href="javascript:view_annotations('+"'"+row["Model"]+"'"+')">'+row["Model"]+"</a>"
        json_str = '{"data":'+table.to_json(orient='records')+'}'
        with open(self.working_dir+"/html/data.json", 'w') as f:
            f.write(json_str)