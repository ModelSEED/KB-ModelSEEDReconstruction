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
from modelseedpy import AnnotationOntology, MSPackageManager,MSGenome, MSMedia, MSModelUtil, MSBuilder, MSGapfill, FBAHelper, MSGrowthPhenotypes, MSModelUtil, MSATPCorrection,MSModelReport
from modelseedpy.helpers import get_template
from modelseedpy.core.msgenomeclassifier import MSGenomeClassifier
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msgenome import normalize_role
from kbbasemodules.basemodelingmodule import BaseModelingModule

logger = logging.getLogger(__name__)

class ModelSEEDRecon(BaseModelingModule):
    def __init__(self,config,module_dir="/kb/module",working_dir=None,token=None,clients={},callback=None):
        BaseModelingModule.__init__(self,"ModelSEEDReconstruction",config,module_dir,working_dir,token,clients,callback)
        self.core_template = None
        self.gs_template = None
        self.version = "0.1.1.msr"
        self.module_dir = module_dir
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        
    def build_metabolic_models(self,params):
        default_media = "KBaseMedia/AuxoMedia"
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
            "gapfilling_media_list":None,
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
            "save_report_to_kbase":True,
            "change_to_complete":False,
            "gapfilling_mode":"Cumulative",
            "base_media":None,
            "compound_list":None,
            "base_media_target_element":"C"
        })
        if params["change_to_complete"]:
            default_media = "KBaseMedia/Complete"
        params["genome_refs"] = self.process_genome_list(params["genome_refs"],params["workspace"])
        #Processing media
        params["gapfilling_media_objs"] = self.process_media_list(params["gapfilling_media_list"],default_media,params["workspace"])
        #Preloading core and preselected template
        self.gs_template = None
        if params["gs_template_ref"]:
            self.gs_template = self.get_template(params["gs_template_ref"],None)
        if params["core_template_ref"]:
            self.core_template = self.get_template(params["core_template_ref"],None)
        else:
            self.core_template = self.get_template(self.templates["core"],None)  
        #Initializing classifier
        genome_classifier = self.get_classifier()
        #Initializing output data tables
        result_table = pd.DataFrame({})
        default_output = {"Model":None,"Genome":None,"Genes":None,"Class":None,
                          "Model genes":None,"Reactions":None,
                          "Core GF":None,"GS GF":None,"Growth":None,"Comments":[]}
        #Retrieving genomes and building models one by one
        mdllist = []
        for i,gen_ref in enumerate(params["genome_refs"]):
            template_type = params["gs_template"]
            #Getting RAST annotated genome, which will be reannotated as needed
            genome = self.get_msgenome_from_ontology(gen_ref,native_python_api=True,output_ws=params["workspace"])
            #Initializing output row
            current_output = default_output.copy()
            current_output["Comments"] = []
            gid = genome.id
            current_output["Model"] = gid+params["suffix"]+'<br><a href="'+gid+params["suffix"]+'-recon.html" target="_blank">(see reconstruction report)</a><br><a href="'+gid+params["suffix"]+'-full.html" target="_blank">(see full view)</a>'
            current_output["Genome"] = genome.info[10]["Name"]
            current_output["Genes"] = genome.info[10]["Number of Protein Encoding Genes"]
            print(genome)
            #Pulling annotation priority
            current_output["Comments"].append("Other annotation priorities not supported by this app yet. Using RAST.")
            if template_type == "auto":
                current_output["Class"] = genome_classifier.classify(genome)
                print(current_output["Class"])
                if current_output["Class"] == "P":
                    current_output["Class"] = "Gram Positive"
                    template_type = "gp"
                elif current_output["Class"] == "N" or current_output["Class"] == "--":
                    current_output["Class"] = "Gram Negative"
                    template_type = "gn"
                elif current_output["Class"] == "A":
                    current_output["Class"] = "Archaea"
                    template_type = "ar"
                elif current_output["Class"] == "C":
                    current_output["Class"] = "Cyanobacteria"
                    template_type = "cyano"
                    current_output["Comments"].append("Cyanobacteria not yet supported. Skipping genome.")
                    result_table = result_table.append(current_output, ignore_index = True)
                    continue
                else:
                    current_output["Comments"].append("Unrecognized genome class "+current_output["Class"]+". Skipping genome.")
                    result_table = result_table.append(current_output, ignore_index = True)
                    continue
            if not self.gs_template:
                self.gs_template = self.get_template(self.templates[template_type],None)
            #Building model            
            base_model = FBAModel({'id':gid+params["suffix"], 'name':genome.scientific_name})
            builder = MSBuilder(genome, self.gs_template)
            mdl = builder.build(base_model, '0', False, False)
            mdl.genome = genome
            mdl.template = self.gs_template
            mdl.core_template_ref = str(self.core_template.info)
            mdl.genome_ref = self.wsinfo_to_ref(genome.info)
            mdl.template_ref = str(self.gs_template.info)
            current_output["Core GF"] = "NA" 
            mdlutl = MSModelUtil.get(mdl)
            if params["atp_safe"]:
                atpcorrection = MSATPCorrection(mdlutl,self.core_template,params["atp_medias"],load_default_medias=params["load_default_medias"],max_gapfilling=params["max_gapfilling"],gapfilling_delta=params["gapfilling_delta"],forced_media=params["forced_atp_list"],default_media_path=self.module_dir+"/data/atp_medias.tsv")
                tests = atpcorrection.run_atp_correction()
                current_output["Core GF"] = len(atpcorrection.cumulative_core_gapfilling)
            #Setting the model ID so the model is saved with the correct name in KBase
            mdlutl.get_attributes()["class"] = current_output["Class"]
            mdlutl.wsid = gid+params["suffix"]
            #Running gapfilling
            current_output["GS GF"] = "NA"
            if params["run_gapfilling"]:
                self.gapfill_metabolic_models({
                    "media_objs":params["gapfilling_media_objs"],#
                    "model_objs":[mdlutl],#
                    "atp_safe":params["atp_safe"],#
                    "workspace":params["workspace"],#
                    "suffix":"",#
                    "default_objective":"bio1",#
                    "output_data":{mdlutl:current_output},#
                    "forced_atp_list":params["forced_atp_list"],
                    "templates":[self.gs_template],
                    "internal_call":True,
                    "gapfilling_mode":params["gapfilling_mode"],
                    "base_media":params["base_media"],
                    "compound_list":params["compound_list"],
                    "base_media_target_element":params["base_media_target_element"]
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
        self.build_dataframe_report(result_table,mdllist)
        if params["save_report_to_kbase"]:
            output = self.save_report_to_kbase()
        if params["return_data"]:
            output["data"] = result_table.to_json()
        if params["return_model_objects"]:
            output["model_objs"] = mdllist
        return output
    
    def gapfill_metabolic_models(self,params):
        default_media = "KBaseMedia/AuxoMedia"
        self.initialize_call("gapfill_metabolic_models",params,True)
        self.validate_args(params,["workspace"],{
            "media_list":None,
            "media_objs":None,
            "model_list":None,
            "model_objectives":[],
            "model_objs":[],
            "atp_safe":True,
            "suffix":".gf",
            "forced_atp_list":[],
            "templates":None,
            "core_template_ref":None,
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
            "return_model_objects":False,
            "return_data":False,
            "save_report_to_kbase":True,
            "change_to_complete":False,
            "gapfilling_mode":"Cumulative",
            "base_media":None,
            "compound_list":None,
            "base_media_target_element":"C"
        })
        base_comments = []
        if params["change_to_complete"]:
            base_comments.append("Changing default to complete.")
            default_media = "KBaseMedia/Complete"
        result_table = pd.DataFrame({})
        default_output = {"Model":None,"Genome":None,"Genes":None,"Class":None,
            "Model genes":None,"Reactions":None,"Core GF":None,"GS GF":None,"Growth":None,"Comments":[]}
        #Retrieving models if not provided already
        if "model_objs" not in params or len(params["model_objs"]) == 0:
            params["model_objs"] = []
            for mdl_ref in params["model_list"]:
                params["model_objs"].append(self.get_model(mdl_ref))
        #Processing media
        if not params["media_objs"]:
            params["media_objs"] = self.process_media_list(params["media_list"],default_media,params["workspace"])
        #Processing compound list
        if params["compound_list"]:
            if not params["base_media"]:
                base_comments.append("No base media provided. Ignoring compound list.")
            else:
                for cpd in params["compound_list"]:
                    newmedia = MSMedia.from_dict({cpd:100})
                    newmedia.merge(params["base_media"])
                    params["media_objs"].append(newmedia)
        #Compiling additional tests
        additional_tests = []
        for i,limit_media in enumerate(params["limit_medias"]):
            additional_tests.append({
                "objective":params["limit_objectives"][i],
                "media":self.get_media(limit_media,None),
                "is_max_threshold":params["is_max_limits"][i],
                "threshold":params["limit_thresholds"][i]
            })
        #Getting core template
        if params["core_template_ref"]:
            self.core_template = self.get_template(params["core_template_ref"],None)
        else:
            self.core_template = self.get_template(self.templates["core"],None)
        #Iterating over each model and running gapfilling
        for i,mdlutl in enumerate(params["model_objs"]):
            current_output = default_output.copy()
            current_output["Comments"] = base_comments
            current_output["Model"] = mdlutl.wsid+params["suffix"]+'<br><a href="'+mdlutl.wsid+params["suffix"]+'-recon.html" target="_blank">(see reconstruction report)</a><br><a href="'+mdlutl.wsid+params["suffix"]+'-full.html" target="_blank">(see full view)</a>'
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
                tests = mdlutl.get_atp_tests(core_template=self.core_template,atp_media_filename=self.module_dir+"/data/atp_medias.tsv",recompute=False)
                additional_tests.extend(tests)
            #Creating gapfilling object and configuring solver
            #mdlutl.model.solver = config["solver"]
            if not params["templates"]:
                params["templates"] = [self.get_template(mdlutl.model.template_ref)]
            msgapfill = MSGapfill(
                mdlutl,
                params["templates"],
                params["source_models"],
                additional_tests,
                blacklist=params["reaction_exlusion_list"],
                default_target=params["model_objectives"][i],
                minimum_obj=params["minimum_objective"],
                base_media=params["base_media"],
                base_media_target_element=params["base_media_target_element"]
            )
            #Running gapfilling in all conditions
            mdlutl.gfutl.cumulative_gapfilling = []
            growth_array = []
            solutions = msgapfill.run_multi_gapfill(
                params["media_objs"],
                target=params["model_objectives"][i],
                default_minimum_objective=params["minimum_objective"],
                binary_check=False,
                prefilter=True,
                check_for_growth=True,
                gapfilling_mode=params["gapfilling_mode"],
                run_sensitivity_analysis=True,
                integrate_solutions=True
            )
            for media in params["media_objs"]:
                if media.id in solutions and "growth" in solutions[media]:
                    growth_array.append(media.id+":"+str(solutions[media]["growth"]))
            current_output["Growth"] = "<br>".join(growth_array)
            current_output["GS GF"] = len(mdlutl.gfutl.cumulative_gapfilling)
            current_output["Reactions"] = mdlutl.nonexchange_reaction_count()
            current_output["Model genes"] = len(mdlutl.model.genes)
            #Saving completely gapfilled model
            self.save_model(mdlutl,params["workspace"],None,params["suffix"])
            if not params["internal_call"]:
                result_table = result_table.append(current_output, ignore_index = True)
        output = {}
        if not params["internal_call"]:
            self.build_dataframe_report(result_table,params["model_objs"])
            if params["save_report_to_kbase"]:
                output = self.save_report_to_kbase()
            if params["return_data"]:
                output["data"] = result_table.to_json()
            if params["return_model_objects"]:
                output["model_objs"] = params["model_objs"]
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
            for model in model_objs:
                msmodrep = MSModelReport(model)
                msmodrep.build_report(self.working_dir+"/html/"+model.wsid+"-recon.html")
                msmodrep.build_multitab_report(self.working_dir+"/html/"+model.wsid+"-full.html") 
        with open(self.working_dir+"/html/index.html", 'w') as f:
            f.write(html)
        #Creating data table file
        json_str = '{"data":'+table.to_json(orient='records')+'}'
        with open(self.working_dir+"/html/data.json", 'w') as f:
            f.write(json_str)