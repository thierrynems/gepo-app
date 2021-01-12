from django.http import HttpResponse, Http404,JsonResponse
from django.shortcuts import render
from datetime import datetime
from .forms import PredictionForm
import requests
from django.core import serializers
import pandas
import sys
import os
from werkzeug.utils import secure_filename
#variable globales 
allSpecies=dict()
allProtein=dict()
#content list of model
urlModel="http://3.133.132.243:9082/deeply/getModels"
MODEL_DICT=dict()
try: 
    r = requests.get(urlModel)
    listModel=r.json()
    for key, oneModel in listModel.items():
        isdisplay=oneModel['isdisplay']
        if(int(isdisplay)==1): 
            model_code=str(oneModel['model_code'])
            model_name=str(oneModel['model_name'])
            MODEL_DICT[model_code]=model_name
except: 
    MODEL_DICT=dict()

#modelDict={"deeplyTrain9deeply-model":"Deeply Core","sgsTrain9deeply-model":"Deeply SGSv0.1","deepHE9deeply-model":"DeeplyHE","SRBTRain9deeply-model":"Deeply Core SRB","SRBSgsTrain9deeply-model":"SGS SRB","SRBHETrain9deeply-model":"DeeplyHE SRB","mlModel":"Classical Machine Learning Model"}

UPLOAD_FOLDER = 'static/model'
ALLOWED_EXTENSIONS = {'h5', 'H5', 'hdf5', 'HDF5'}
#URL des API
URL_API="http://3.133.132.243:5001"
#list of species
URL_API_SPECIES=str(URL_API+"/searchspecies")
URL_API_SPECIES_KEGG=str(URL_API+"/searchkeggorganis")
#protein 
URL_API_PROTEIN_EPATH=str(URL_API+"/searchproteinepath")

URL_API_PROTEIN_STRING=str(URL_API+"/searchproteinbyname")




def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

#vue de la page d'accueil
def index(request):
    return render(request, 'blog/index.html', locals())

#vue du formulaire de prediction
def predictform(request):
    modelDict=MODEL_DICT
    #modelDict={"deeply":"Deeply Essential Model","mlModel":"Classical Machine Learning Model"}
    #modelDict={"deeplyTrain9deeply-model":"Deeply Core","sgsTrain9deeply-model":"Deeply SGSv0.1","deepHE9deeply-model":"DeeplyHE","SRBTRain9deeply-model":"Deeply Core SRB","SRBSgsTrain9deeply-model":"SGS SRB","SRBHETrain9deeply-model":"DeeplyHE SRB","mlModel":"Classical Machine Learning Model"}
    actionDict={"feature":"Feature Engeneering","prediction":"Prediction"}
    defaultGeneSeq="ATGGAAAATATATTAGACCTGTGGAACCAAGCCCTTGCTCAAATCGAAAAAAAGTTGAGCAAACCGAGTTTTGAGACTTGGATGAAGTCAACCAAAGCCCACTCACTGCAAGGCGATACATTAACAATCACGGCTCCCAATGAATTTGCCAGAGACTGGCTGGAGTCCAGATACTTGCATCTGATTGCAGATACTATATATGAATTAACCGGGGAAGAATTGAGCATTAAGTTTGTCATTCCTCAAAATCAAGATGTTGAGGACTTTATGCCGAAACCGCAAGTCAAAAAAGCGGTCAAAGAAGATACATCTGATTTTCCTCAAAATATGCTCAATCCAAAATATACTTTTGATACTTTTGTCATCGGATCTGGAAACCGATTTGCACATGCTGCTTCCCTCGCAGTAGCGGAAGCGCCCGCGAAAGCTTACAACCCTTTATTTATCTATGGGGGCGTCGGCTTAGGGAAAACACACTTAATGCATGCGATCGGCCATTATGTAATAGATCATAATCCTTCTGCCAAAGTGGTTTATCTGTCTTCTGAGAAATTTACAAACGAATTCATCAACTCTATCCGAGATAATAAAGCCGTCGACTTCCGCAATCGCTATCGAAATGTTGATGTGCTTTTGATAGATGATATTCAATTTTTAGCGGGGAAAGAACAAACCCAGGAAGAATTTTTCCATACATTTAACACATTACACGAAGAAAGCAAACAAATCGTCATTTCAAGTGACCGGCCGCCAAAGGAAATTCCGACACTTGAAGACAGATTGCGCTCACGTTTTGAATGGGGACTTATTACAGATATCACACCGCCTGATCTAGAAACGAGAATTGCAATTTTAAGAAAAAAGGCCAAAGCAGAGGGCCTCGATATTCCGAACGAGGTTATGCTTTACATCGCGAATCAAATCGACAGCAATATTCGGGAACTCGAAGGAGCATTAATCAGAGTTGTCGCTTATTCATCTTTAATTAATAAAGATATTAATGCTGATCTGGCCGCTGAGGCGTTGAAAGATATTATTCCTTCCTCAAAACCGAAAGTCATTACGATAAAAGAAATTCAGAGGGTAGTAGGCCAGCAATTTAATATTAAACTCGAGGATTTCAAAGCAAAAAAACGGACAAAGTCAGTAGCTTTTCCGCGTCAAATCGCCATGTACTTATCAAGGGAAATGACTGATTCCTCTCTTCCTAAAATCGGTGAAGAGTTTGGAGGACGTGATCATACGACCGTTATTCATGCGCATGAAAAAATTTCAAAACTGCTGGCAGATGATGAACAGCTTCAGCAGCATGTAAAAGAAATTAAAGAACAGCTTAAATAG"
    defaultProtSeq="MENILDLWNQALAQIEKKLSKPSFETWMKSTKAHSLQGDTLTITAPNEFARDWLESRYLHLIADTIYELTGEELSIKFVIPQNQDVEDFMPKPQVKKAVKEDTSDFPQNMLNPKYTFDTFVIGSGNRFAHAASLAVAEAPAKAYNPLFIYGGVGLGKTHLMHAIGHYVIDHNPSAKVVYLSSEKFTNEFINSIRDNKAVDFRNRYRNVDVLLIDDIQFLAGKEQTQEEFFHTFNTLHEESKQIVISSDRPPKEIPTLEDRLRSRFEWGLITDITPPDLETRIAILRKKAKAEGLDIPNEVMLYIANQIDSNIRELEGALIRVVAYSSLINKDINADLAAEALKDIIPSSKPKVITIKEIQRVVGQQFNIKLEDFKAKKRTKSVAFPRQIAMYLSREMTDSSLPKIGEEFGGRDHTTVIHAHEKISKLLADDEQLQQHVKEIKEQLK"
    #text="bonjour"

    # Quoiqu'il arrive, on affiche la page du formulaire.
    return render(request, 'blog/predictForm.html', locals())

#soumission ajax
def ajaxpredictform(request):
    
    if request.method == 'POST':
        modelSelect = request.POST.get('modelSelect', None)
        actionSelect = request.POST.get('actionSelect', None)
        geneSequence = request.POST.get('geneSequence', None)
        proteinSequence = request.POST.get('proteinSequence', None)
        feature=False
        prediction=False
        envoi=True
        modelDict=MODEL_DICT
        #modelDict={"deeply":"Deeply Essential Model","mlModel":"Classical Machine Learning Model"}

        #modelDict={"deeplyTrain9deeply-model":"Deeply Core","sgsTrain9deeply-model":"Deeply SGSv0.1","deepHE9deeply-model":"DeeplyHE","SRBTRain9deeply-model":"Deeply Core SRB","SRBSgsTrain9deeply-model":"SGS SRB","SRBHETrain9deeply-model":"DeeplyHE SRB","mlModel":"Classical Machine Learning Model"}
        actionDict={"feature":"Feature Engeneering","prediction":"Prediction"}
        data = {'geneSeq': geneSequence, 'protSeq': proteinSequence,'model':modelSelect}
        args = {'geneSeq': geneSequence, 'protSeq': proteinSequence,'model':modelSelect}
        modelChoice=modelDict[modelSelect]
        header = {'Content-type': 'application/x-www-form-urlencoded; charset=utf-8'}
        service="" #le service choisit (feature engeneering or prediction)
        isEssential=False # test si le gene est esentiel
        inConstruct=False #modele en construction 
        if(modelSelect=="mlModel"):inConstruct=True
        if(actionSelect=="feature"):
            feature=True
            service=actionDict['feature']
            #url="http://ec2-18-222-229-153.us-east-2.compute.amazonaws.com/deeply/feature"
            url="http://3.133.132.243:9082/deeply/feature"
            r = requests.post(url, data=data, params=args, headers=header, verify=False)
            feauturePredict=r.json()
            CIARCSUFeat=feauturePredict['CIARCSUFeat']['DEG100011']
            GCFeat=feauturePredict['GCFeat']['DEG100011']
            GeneLength=feauturePredict['GeneLength']['DEG100011']
            KmerFeat=feauturePredict['KmerFeat']['DEG100011']
            ProteinFeat=feauturePredict['ProteinFeat']['DEG100011']
                
        else:
            prediction=True
            service=actionDict['prediction']
            #url="http://ec2-18-222-229-153.us-east-2.compute.amazonaws.com/deeply/prediction"
            url="http://3.133.132.243:9082/deeply/prediction"
            r = requests.post(url, data=data, params=args, headers=header, verify=False)
            resultPredict=r.json()
            EssentialGeneProbability=resultPredict['EssentialGeneProbability']
            NonEssentialGeneProbability=resultPredict['NonEssentialGeneProbability']
            essentialProb=EssentialGeneProbability.replace("[",'')
            essentialProb=essentialProb.replace("]",'')
            if(float(essentialProb)>0.5): isEssential=True
         
        #htmlRender=serializers.serialize("json", htmlRender)
       # response = {'html':htmlRender # response message}

        return render(request, 'blog/ajaxpredictForm.html', locals())
            
#vue du formulaire recherche par protein
#seach by name 
def searchbyprotein(request):
    #modelDict={"deeply":"Deeply Essential Model","mlModel":"Classical Machine Learning Model"}
    modelDict=MODEL_DICT
    #modelDict={"deeplyTrain9deeply-model":"Deeply Core","sgsTrain9deeply-model":"Deeply SGSv0.1","deepHE9deeply-model":"DeeplyHE","SRBTRain9deeply-model":"Deeply Core SRB","SRBSgsTrain9deeply-model":"SGS SRB","SRBHETrain9deeply-model":"DeeplyHE SRB","mlModel":"Classical Machine Learning Model"}
    actionDict={"feature":"Feature Engeneering","prediction":"Prediction"}
    defaultGeneSeq="ATGGAAAATATATTAGACCTGTGGAACCAAGCCCTTGCTCAAATCGAAAAAAAGTTGAGCAAACCGAGTTTTGAGACTTGGATGAAGTCAACCAAAGCCCACTCACTGCAAGGCGATACATTAACAATCACGGCTCCCAATGAATTTGCCAGAGACTGGCTGGAGTCCAGATACTTGCATCTGATTGCAGATACTATATATGAATTAACCGGGGAAGAATTGAGCATTAAGTTTGTCATTCCTCAAAATCAAGATGTTGAGGACTTTATGCCGAAACCGCAAGTCAAAAAAGCGGTCAAAGAAGATACATCTGATTTTCCTCAAAATATGCTCAATCCAAAATATACTTTTGATACTTTTGTCATCGGATCTGGAAACCGATTTGCACATGCTGCTTCCCTCGCAGTAGCGGAAGCGCCCGCGAAAGCTTACAACCCTTTATTTATCTATGGGGGCGTCGGCTTAGGGAAAACACACTTAATGCATGCGATCGGCCATTATGTAATAGATCATAATCCTTCTGCCAAAGTGGTTTATCTGTCTTCTGAGAAATTTACAAACGAATTCATCAACTCTATCCGAGATAATAAAGCCGTCGACTTCCGCAATCGCTATCGAAATGTTGATGTGCTTTTGATAGATGATATTCAATTTTTAGCGGGGAAAGAACAAACCCAGGAAGAATTTTTCCATACATTTAACACATTACACGAAGAAAGCAAACAAATCGTCATTTCAAGTGACCGGCCGCCAAAGGAAATTCCGACACTTGAAGACAGATTGCGCTCACGTTTTGAATGGGGACTTATTACAGATATCACACCGCCTGATCTAGAAACGAGAATTGCAATTTTAAGAAAAAAGGCCAAAGCAGAGGGCCTCGATATTCCGAACGAGGTTATGCTTTACATCGCGAATCAAATCGACAGCAATATTCGGGAACTCGAAGGAGCATTAATCAGAGTTGTCGCTTATTCATCTTTAATTAATAAAGATATTAATGCTGATCTGGCCGCTGAGGCGTTGAAAGATATTATTCCTTCCTCAAAACCGAAAGTCATTACGATAAAAGAAATTCAGAGGGTAGTAGGCCAGCAATTTAATATTAAACTCGAGGATTTCAAAGCAAAAAAACGGACAAAGTCAGTAGCTTTTCCGCGTCAAATCGCCATGTACTTATCAAGGGAAATGACTGATTCCTCTCTTCCTAAAATCGGTGAAGAGTTTGGAGGACGTGATCATACGACCGTTATTCATGCGCATGAAAAAATTTCAAAACTGCTGGCAGATGATGAACAGCTTCAGCAGCATGTAAAAGAAATTAAAGAACAGCTTAAATAG"
    defaultProtSeq="MENILDLWNQALAQIEKKLSKPSFETWMKSTKAHSLQGDTLTITAPNEFARDWLESRYLHLIADTIYELTGEELSIKFVIPQNQDVEDFMPKPQVKKAVKEDTSDFPQNMLNPKYTFDTFVIGSGNRFAHAASLAVAEAPAKAYNPLFIYGGVGLGKTHLMHAIGHYVIDHNPSAKVVYLSSEKFTNEFINSIRDNKAVDFRNRYRNVDVLLIDDIQFLAGKEQTQEEFFHTFNTLHEESKQIVISSDRPPKEIPTLEDRLRSRFEWGLITDITPPDLETRIAILRKKAKAEGLDIPNEVMLYIANQIDSNIRELEGALIRVVAYSSLINKDINADLAAEALKDIIPSSKPKVITIKEIQRVVGQQFNIKLEDFKAKKRTKSVAFPRQIAMYLSREMTDSSLPKIGEEFGGRDHTTVIHAHEKISKLLADDEQLQQHVKEIKEQLK"
    #url="http://195.24.221.74:9200/species*/_search?size=150&q=*.*&pretty=true"
    #api to get list of species
    urlAWS="http://3.133.132.243:9082/deeply/speciesList"

    #url="http://195.24.221.74:5001/searchspecies"
    #r = requests.post(url, data=data, params=args, headers=header, verify=False)
    try:
        r = requests.get(URL_API_SPECIES)
        species=r.json()
        allSpecies=species
    except:
        try: 
            r = requests.get(urlAWS)
            species=r.json()
            allSpecies=species
        except:
            species=dict()

    #load dataset 
    #dataset="blog/datasets/species.v11.0.txt"
    #df = pandas.read_csv(dataset,sep= '\t', header = 0)

    #row and col
    #row,col = df.shape
    #init dictionary of organism
    organismDict=dict()
    #mySpecies=species['hits']['hits']
    #rows =len(mySpecies)
    
    num=0

    for key, oneSpecies in species.items():
    	try:
    		specieId=oneSpecies['species_id']
    		specieName=oneSpecies['official_name']
    		organismDict[specieId]=specieName
    	except Exception as  e: 
    		print(e)


    #for index in range(row) :
    #    taxaId=df['taxon_id'][index]
    #    organism=df['official_name_NCBI'][index]
    #    organismDict[taxaId]=organism

    #text="bonjour"

    # Quoiqu'il arrive, on affiche la page du formulaire.
    return render(request, 'blog/searchbyprotein.html', locals())

#ajax submit protien search by name 
#soumission ajax
def ajaxsearchprotein(request):
    
    if request.method == 'POST':
        action = request.POST.get('action', None)
        #file to download 
        baseDir="static/download/"
        #test if baseDir exists
        if os.path.isdir(baseDir): 
            print("le repertoire existe")
        else: 
            print("Creation du repertoire "+ baseDir)
            os.mkdir(baseDir);
        if action =='step1': 
            proteinName = request.POST.get('proteinName', None)
            specie = request.POST.get('specie', None)
            step=action
            dataaws = {'proteinName': proteinName,'species_id':specie}
            argsaws = {'proteinName': proteinName,'species_id':specie}
            data = {'protein_name': proteinName,'species_id':specie}
            args = {'protein_name': proteinName,'species_id':specie}
            actionDict={"feature":"Feature Engeneering","prediction":"Prediction"}
            header = {'Content-type': 'application/x-www-form-urlencoded; charset=utf-8'}
            urlAWS="http://3.133.132.243:9082/deeply/searchProtein"
            url="http://195.24.221.74:5001/searchproteinbyname"
            #request to get list of species share the protein name
            proteinList=dict()
            try:
                r = requests.get(URL_API_PROTEIN_STRING, params=args) #urlAWS
                proteinList=r.json()
            except:
                try: 
                    r = requests.post(urlAWS, data=dataaws, params=argsaws, headers=header, verify=False) #url
                    proteinList=r.json()
                except Exception as e:
                    proteinList=dict()

            #get protein ID
            proteinID=""
            proteinExternalId=""
            for key, protein in proteinList.items(): 
                proteinID=protein['protein_id']
                proteinExternalId=protein['protein_external_id']
                break

            filename="resultatStep1"
            currentfileCSV=os.path.join(baseDir, filename+'.csv')
            currentfileTXT=os.path.join(baseDir, filename+'.txt')
            headercsv="#,TaxaId,Species,Annotation,Protein Name,Protein Size,Protein External ID\n"
            headertxt="#    taxaId  Species Annotation  Protein Name  Protein Size  Protein External ID\n"
            filecsv = open(currentfileCSV, "w")
            filetxt = open(currentfileTXT, "w")
            filecsv.write(str(headercsv))
            filetxt.write(str(headertxt))
            lineIndex=1
            for key, protein in proteinList.items(): 
                proteinID=str(protein['protein_id'])
                taxaId=str(protein['species_id'])
                protein_size=str(protein['protein_size'])
                annotation=str(protein['annotation'])
                specie_name=str(protein['species_compact_name'])
                #specie_name=str(protein['specie_name'])
                proteinExternalId=str(protein['protein_external_id'])
                protein_name=str(protein['preferred_name'])
                linecsv=str(lineIndex)+","+taxaId+","+specie_name+","+annotation+","+protein_name+","+protein_size+","+proteinExternalId+"\n"
                linetxt=str(lineIndex)+"    "+taxaId+"  "+specie_name+" "+annotation+"  "+protein_name+"    "+protein_size+"    "+proteinExternalId+"\n"
                filecsv.write(linecsv)
                filetxt.write(linetxt)
                lineIndex=lineIndex+1
            filecsv.close()
            filetxt.close()
            # resquest to get list of organism
            #url="http://3.133.132.243:9082/deeply/speciesList"
            #r = requests.get(url)
            #species=r.json()
            return render(request, 'blog/ajaxsearchprotein.html', locals())
        elif action == 'step2': 
            proteinName = request.POST.get('proteinName', None)
            specie = request.POST.get('specie', None)
            service= request.POST.get('service', None)
            listCheck= request.POST.get('listCheck', None)
            proteinID= request.POST.get('proteinID', None)
            proteinExternalId=request.POST.get('proteinExternalId', None)
            modelDict=MODEL_DICT
            #modelDict={"deeplyTrain9deeply-model":"Deeply Core","sgsTrain9deeply-model":"Deeply SGSv0.1","deepHE9deeply-model":"DeeplyHE","SRBTRain9deeply-model":"Deeply Core SRB","SRBSgsTrain9deeply-model":"SGS SRB","SRBHETrain9deeply-model":"DeeplyHE SRB"}
            #modelDict={"deeplyTrain9deeply-model":"DeeplyCore","sgsTrain9deeply-model":"SGSv0.1","deepHE9deeply-model":"DeeplyHE","SRBTRain9deeply-model":"SRB-DeeplyCore","SRBSgsTrain9deeply-model":"SRB-SGSv0.1","SRBHETrain9deeply-model":"SRB-DeeplyHE"}
            #"mlModel":"Classical Machine Learning Model"
            actionDict={"feature":"Feature Engeneering","prediction":"Prediction"}

            #obtain protein_external_id : structure: 756067.MicvaDRAFT_3438
            #obtain: list of species check 
            urlSpeciesKegg="http://3.133.132.243:9082/deeply/searchKeggOrganism"
            header = {'Content-type': 'application/x-www-form-urlencoded; charset=utf-8'}
            listSpecies=listCheck.split(sep=";")
            index=0
            resultDict=dict()
            urlSpecies="http://3.133.132.243:9082/deeply/speciesList"
            urlFeature="http://3.133.132.243:9082/deeply/feature"
            urlPredict="http://3.133.132.243:9082/deeply/prediction"
            #r = requests.post(url, data=data, params=args, headers=header, verify=False)
            #r = requests.get(urlSpecies) #URL_API_SPECIES
            r = requests.get(URL_API_SPECIES) #URL_API_SPECIES
            speciesAll=r.json()     
            numberLocus=""
            isNumber=False        
            #resultDict[index]=item
            #index=index+1

            filenamePredict="predictResult"
            filenameFeature="FeatureResult"
            ########################################file prediction 
            currentfilePredictCSV=os.path.join(baseDir, filenamePredict+'.csv')
            currentfilePredictTXT=os.path.join(baseDir, filenamePredict+'.txt')
            headercsv="#,TaxaId,Species,gene_locus,Protein Name,Model,E_score,NE_score,Decision,NTseq,AAseq\n"
            headertxt="#    TaxaId    Species    gene_locus    Protein Name    Model    E_score    NE_score    Decision    NTseq    AAseq\n"
            filePredictcsv = open(currentfilePredictCSV, "w")
            filePredicttxt = open(currentfilePredictTXT, "w")
            filePredictcsv.write(str(headercsv))
            filePredicttxt.write(str(headertxt))
            ################################################feature
            currentfileFeatCSV=os.path.join(baseDir, filenameFeature+'.csv')
            currentfileFeatTXT=os.path.join(baseDir, filenameFeature+'.txt')
            fileFeatcsv = open(currentfileFeatCSV, "w")
            fileFeattxt = open(currentfileFeatTXT, "w")
            columnFeat=['NT-Length','GC','CIA','RCSU','TTT', 'TTC', 'TTA', 'TTG', 'CTT','CTC', 'CTA', 'CTG', 'ATT', 'ATC','ATA', 'ATG', 'GTT', 'GTC', 'GTA','GTG', 'TAT', 'TAC', 'TAA', 'TAG','CAT', 'CAC', 'CAA', 'CAG', 'AAT','AAC', 'AAA', 'AAG', 'GAT', 'GAC','GAA', 'GAG', 'TCT', 'TCC', 'TCA','TCG', 'CCT', 'CCC', 'CCA', 'CCG','ACT', 'ACC', 'ACA', 'ACG', 'GCT','GCC', 'GCA', 'GCG', 'TGT', 'TGC','TGA', 'TGG', 'CGT', 'CGC', 'CGA','CGG', 'AGT', 'AGC', 'AGA', 'AGG','GGT', 'GGC', 'GGA', 'GGG','A', 'R', 'N', 'D','C', 'Q', 'E', 'G','H', 'I', 'L', 'K','M', 'F', 'P', 'S','T', 'W', 'Y', 'V','AA-Length'];
            headerFeatcsv=",".join(columnFeat)+"\n"
            headerFeattxt=" ".join(columnFeat)+"\n"
            fileFeatcsv.write(str(headerFeatcsv))
            fileFeattxt.write(str(headerFeattxt))

            lineIndex=1
            numberLine=1
            for speciesId in listSpecies : 
                for key, species in speciesAll.items():
                    if  int(species['species_id']) == int(speciesId) :
                        item=dict() 
                        officialName=species['official_name'] 
                        compactName=species['compact_name']
                        is_sequence=True
                        data = {'officialName': officialName,"compactName":compactName,"proteinName":proteinName}
                        args = {'officialName': officialName,"compactName":compactName,"proteinName":proteinName}
                        try: 
                            r = requests.post(urlSpeciesKegg, data=data, params=args, headers=header, verify=False)
                            speciesInfo=r.json()
                            organism_code=speciesInfo['organism_code']
                            organism_taxa=speciesInfo['organism_taxa']
                            geneLocus=speciesInfo['gene_locus']
                        except:
                            organism_code=""
                            organism_taxa=""
                            geneLocus=""
                        
                        ##################################get sequene ################################
                        subUrl=str(organism_code)+":"+str(geneLocus)
                        data = {'subURL': subUrl,'locusNumber':numberLocus,"organism_code":organism_code,"isNumber":isNumber}
                        arg = {'subURL': subUrl,'locusNumber':numberLocus,"organism_code":organism_code,"isNumber":isNumber}
                        urlSequence="http://3.133.132.243:9082/deeply/getSequence"
                        try: 
                            r = requests.post(urlSequence, data=data, params=args, headers=header, verify=False)
                            sequences=r.json()
                            sequenceAA=sequences["aa"]
                            sequenceNT=sequences["nt"]
                        except:
                            sequenceAA=""
                            sequenceNT=""
                        ###########################
                        allModel=dict()
                        featureDict=dict()
                        if sequenceAA=="" and sequenceNT=="":
                            sequenceAA="sequence not avalaible"
                            sequenceNT="sequence not avalaible"
                            is_sequence=False
                        else:
                            ##########################feature engineering ####################################
                            lineFeatcsv=""
                            lineFeattxt=""
                            sequenceNT=sequenceNT.upper()
                            sequenceAA=sequenceAA.upper()
                            dataFeature = {'geneSeq': sequenceNT, 'protSeq': sequenceAA}
                            argsFeature = {'geneSeq': sequenceNT, 'protSeq': sequenceAA}
                            try:
                                r = requests.post(urlFeature, data=dataFeature, params=argsFeature, headers=header, verify=False)
                                feauturePredict=r.json()
                                CIARCSUFeat=feauturePredict['CIARCSUFeat']['DEG100011']
                                GCFeat=feauturePredict['GCFeat']['DEG100011']
                                GCFeat=GCFeat[0]
                                GeneLength=feauturePredict['GeneLength']['DEG100011']
                                KmerFeat=feauturePredict['KmerFeat']['DEG100011']
                                ProteinFeat=feauturePredict['ProteinFeat']['DEG100011']
                                lineFeatcsv=str(GeneLength)+","+str(GCFeat)+","+",".join([str(_) for _ in CIARCSUFeat])+","+",".join([str(_) for _ in KmerFeat])+","+",".join([str(_) for _ in ProteinFeat])+"\n"
                                lineFeattxt=str(GeneLength)+"  "+str(GCFeat)+" "+" ".join([str(_) for _ in CIARCSUFeat])+" "+" ".join([str(_) for _ in KmerFeat])+"    "+" ".join([str(_) for _ in ProteinFeat])+"\n"
                                fileFeatcsv.write(str(lineFeatcsv))
                                fileFeattxt.write(str(lineFeattxt))

                            except: 
                                CIARCSUFeat=""
                                GCFeat=""
                                GeneLength=""
                                KmerFeat=""
                                ProteinFeat=""
                                lineFeatcsv=""
                                lineFeattxt=""

                                                        #write to file 
                            #fileFeatcsv.write(str(lineFeatcsv))
                            #fileFeattxt.write(str(lineFeattxt))
                            ###############################predict ##########################################
                            for  modelSelect, valueModel in modelDict.items():
                                modelResult=dict()
                                linetxt=""
                                linecsv=""
                                data = {'geneSeq': sequenceNT, 'protSeq': sequenceAA,'model':modelSelect}
                                args = {'geneSeq': sequenceNT, 'protSeq': sequenceAA,'model':modelSelect}

                                try: 
                                    r = requests.post(urlPredict, data=data, params=args, headers=header, verify=False)
                                    resultPredict=r.json()
                                    EssentialGeneProbability=resultPredict['EssentialGeneProbability']
                                    NonEssentialGeneProbability=resultPredict['NonEssentialGeneProbability']
                                    decision=resultPredict['decision']
                                except: 
                                    EssentialGeneProbability=""
                                    NonEssentialGeneProbability=""
                                    decision=""

                                linetxt=str(numberLine)+"   "+str(organism_taxa)+"   "+str(compactName)+"   "+str(geneLocus)+"  "+str(proteinName)
                                linecsv=str(numberLine)+","+str(organism_taxa)+","+str(compactName)+","+str(geneLocus)+","+str(proteinName)
                                linecsv=linecsv+","+str(valueModel)+","+str(EssentialGeneProbability)+","+str(NonEssentialGeneProbability)+","+str(decision)+","+str(sequenceNT)+","+str(sequenceAA)+"\n"
                                linetxt=linetxt+"   "+str(valueModel)+" "+str(EssentialGeneProbability)+"   "+str(NonEssentialGeneProbability)+"    "+str(decision)+"   "+str(sequenceNT)+" "+str(sequenceAA)+"\n"
                                filePredictcsv.write(str(linecsv))
                                filePredicttxt.write(str(linetxt))
                                modelResult={
                                'CIARCSUFeat':CIARCSUFeat,
                                'GCFeat':GCFeat,
                                'GeneLength':GeneLength,
                                'KmerFeat':KmerFeat,
                                'ProteinFeat':ProteinFeat,
                                'EssentialGeneProbability':EssentialGeneProbability,
                                'NonEssentialGeneProbability':NonEssentialGeneProbability,
                                'decision': decision
                                }
                                featureDict={
                                'CIARCSUFeat':CIARCSUFeat,
                                'GCFeat':GCFeat,
                                'GeneLength':GeneLength,
                                'KmerFeat':KmerFeat,
                                'ProteinFeat':ProteinFeat
                                }
                                allModel[valueModel]=modelResult
                                numberLine=numberLine+1
                        ##############################################################################
                        item={
                        'officialName':officialName,
                        'organism_code':organism_code,
                        'sequenceAA':sequenceAA,
                        'sequenceNT':sequenceNT,
                        'geneLocus': geneLocus,
                        'model':allModel,
                        'modelDict':modelDict,
                        'feature': featureDict,
                        'is_sequence':is_sequence
                        }

                        resultDict[index]=item
                        index=index+1
                        break
            #close file
            fileFeatcsv.close()
            fileFeattxt.close()
            filePredictcsv.close()
            filePredicttxt.close()


            #obtain the sequence for protein
            
            

            return render(request, 'blog/ajaxsearchprotein.html', locals())
#vue du formulaire recherche par plusieurs proteines
def searchlistprotein(request):
    #load dataset 
    url="http://3.133.132.243:9082/deeply/speciesList"
    #r = requests.post(url, data=data, params=args, headers=header, verify=False)
    try: 
        r = requests.get(URL_API_SPECIES)
        species=r.json()
        allSpecies=species 
    except Exception as err:
        try:
            r = requests.get(url)
            species=r.json()
            allSpecies=species 
        except: 
            print(err)

    #load dataset 
    #dataset="blog/datasets/species.v11.0.txt"
    #df = pandas.read_csv(dataset,sep= '\t', header = 0)

    #row and col
    #row,col = df.shape
    #init dictionary of organism
    organismDict=dict()
    for key, oneSpecies in species.items():
        specieId=oneSpecies['species_id']
        specieName=oneSpecies['official_name']
        organismDict[specieId]=specieName
    #text="bonjour"

    # Quoiqu'il arrive, on affiche la page du formulaire.
    return render(request, 'blog/searchlistprotein.html', locals())
#vue ajax search from multiple organism
def ajaxsearchlistprotein(request):
    if request.method == 'POST':
        action = request.POST.get('action', None)
        #file to download 
        baseDir="static/download/"
        #test if baseDir exists
        if os.path.isdir(baseDir): 
            print("le repertoire existe")
        else: 
            print("Creation du repertoire "+ baseDir)
            os.mkdir(baseDir);

        if action =='step1': 
            proteinList = request.POST.get('proteinList', None)
            specie = request.POST.get('organism', None)
            actionDict={"feature":"Feature Engeneering","prediction":"Prediction"}
            step=action
            data = {'proteinList': proteinList,'species_id':specie}
            args = {'proteinList': proteinList,'species_id':specie}
            header = {'Content-type': 'application/x-www-form-urlencoded; charset=utf-8'}
            url="http://3.133.132.243:9082/deeply/getProteinList"
            #request to get list of species share the protein name
            r = requests.post(url, data=data, params=args, headers=header, verify=False)
            listProtein=r.json()
            #number of protein query 
            tabProteinList=proteinList.split(";")
            number_protein_query=len(tabProteinList)
            #write in file 
            filename="resultatListProteinStep1"
            currentfileCSV=os.path.join(baseDir, filename+'resultDict.csv')
            currentfileTXT=os.path.join(baseDir, filename+'.txt')
            headercsv="#,TaxaId,Species,#Protein/Query,Total protein of organism,Mapped item\n"
            headertxt="#    TaxaId    Species    #Protein/Query    Total protein of organism    Mapped item\n"
            filecsv = open(currentfileCSV, "w")
            filetxt = open(currentfileTXT, "w")
            filecsv.write(str(headercsv))
            filetxt.write(str(headertxt))
            lineIndex=1
            for key, protein in listProtein.items():
                taxaId=""
                specie_name=str(protein['specie_name'])
                proteinQuery=str( protein['organism_protein_number'])+"/"+str(number_protein_query)
                totalProtein=""
                mappedItem=str(protein['gene_name'])
                linecsv=str(lineIndex)+","+taxaId+","+specie_name+","+proteinQuery+","+totalProtein+","+mappedItem+"\n"
                linetxt=str(lineIndex)+"    "+taxaId+"  "+specie_name+" "+proteinQuery+"  "+totalProtein+"    "+mappedItem+"\n"
                filecsv.write(linecsv)
                filetxt.write(linetxt)
                lineIndex=lineIndex+1
            filecsv.close()
            filetxt.close()

        elif action =='step2':
            proteinList = request.POST.get('proteinName', None)
            specie = request.POST.get('specie', None)
            service=request.POST.get('service', None)
            listCheck=request.POST.get('listCheck', None)
            #parcours de la liste de gene_locus
            allcheckTab=listCheck.split("|")
            proteinNameTab=allcheckTab[0]
            proteinNameList=proteinNameTab.split(",")
            #rawDataCheck=listCheck.split(":")
            rawDataCheck=allcheckTab[1].split(":")
            organism_code=rawDataCheck[0]
            allgeneLocusRaw=rawDataCheck[1]
            allgeneLocus=allgeneLocusRaw.split(",")
            numberLocus=0
            isNumber=False
            header = {'Content-type': 'application/x-www-form-urlencoded; charset=utf-8'}
            modelDict=MODEL_DICT
            #modelDict={"deeplyTrain9deeply-model":"DeeplyCore","sgsTrain9deeply-model":"SGSv0.1","deepHE9deeply-model":"DeeplyHE","SRBTRain9deeply-model":"SRB-DeeplyCore","SRBSgsTrain9deeply-model":"SRB-SGSv0.1","SRBHETrain9deeply-model":"SRB-DeeplyHE"}
            urlSequence="http://3.133.132.243:9082/deeply/getSequence"
            urlFeature="http://3.133.132.243:9082/deeply/feature"
            urlPredict="http://3.133.132.243:9082/deeply/prediction"
            index=0
            resultDict=dict()
            filenamePredict="predictResultListProtein"
            filenameFeature="FeatureResultListProtein"
            ########################################file prediction 
            currentfilePredictCSV=os.path.join(baseDir, filenamePredict+'.csv')
            currentfilePredictTXT=os.path.join(baseDir, filenamePredict+'.txt')
            headercsv="#,Gene Locus,Protein Name,Model,E_score,NE_score,Decision,NTseq,AAseq\n"
            headertxt="#    Gene Locus    Protein Name    Model    E_score    NE_score    Decision    NTseq    AAseq\n"
            filePredictcsv = open(currentfilePredictCSV, "w")
            filePredicttxt = open(currentfilePredictTXT, "w")
            filePredictcsv.write(str(headercsv))
            filePredicttxt.write(str(headertxt))
            ################################################feature
            currentfileFeatCSV=os.path.join(baseDir, filenameFeature+'.csv')
            currentfileFeatTXT=os.path.join(baseDir, filenameFeature+'.txt')
            fileFeatcsv = open(currentfileFeatCSV, "w")
            fileFeattxt = open(currentfileFeatTXT, "w")
            columnFeat=['NT-Length','GC','CIA','RCSU','TTT', 'TTC', 'TTA', 'TTG', 'CTT','CTC', 'CTA', 'CTG', 'ATT', 'ATC','ATA', 'ATG', 'GTT', 'GTC', 'GTA','GTG', 'TAT', 'TAC', 'TAA', 'TAG','CAT', 'CAC', 'CAA', 'CAG', 'AAT','AAC', 'AAA', 'AAG', 'GAT', 'GAC','GAA', 'GAG', 'TCT', 'TCC', 'TCA','TCG', 'CCT', 'CCC', 'CCA', 'CCG','ACT', 'ACC', 'ACA', 'ACG', 'GCT','GCC', 'GCA', 'GCG', 'TGT', 'TGC','TGA', 'TGG', 'CGT', 'CGC', 'CGA','CGG', 'AGT', 'AGC', 'AGA', 'AGG','GGT', 'GGC', 'GGA', 'GGG','A', 'R', 'N', 'D','C', 'Q', 'E', 'G','H', 'I', 'L', 'K','M', 'F', 'P', 'S','T', 'W', 'Y', 'V','AA-Length'];
            headerFeatcsv=",".join(columnFeat)+"\n"
            headerFeattxt=" ".join(columnFeat)+"\n"
            fileFeatcsv.write(str(headerFeatcsv))
            fileFeattxt.write(str(headerFeattxt))
            numberLine=1
            #position=0
            for gene_locus in allgeneLocus:
                subUrl=str(organism_code)+":"+str(gene_locus)
                errorSeq=""
                errorApi=""
                lineFeatcsv=""
                lineFeattxt=""
                #subUrl=gene_locus
                #tabLocus=gene_locus.split(":")
                #organism_code=tabLocus[0]
                data = {'subURL': subUrl,'locusNumber':numberLocus,"organism_code":organism_code,"isNumber":isNumber}
                args = {'subURL': subUrl,'locusNumber':numberLocus,"organism_code":organism_code,"isNumber":isNumber}
                try: 
                    r = requests.post(urlSequence, data=data, params=args, headers=header, verify=False)
                    sequences=r.json()
                    sequenceAA=sequences["aa"]
                    sequenceNT=sequences["nt"]
                    #convertir en majuscule
                    sequenceNT=sequenceNT.upper()
                    sequenceAA=sequenceAA.upper()
                except Exception as e: 
                    sequenceAA=""
                    sequenceNT=""
                    errorSeq=sequences["error"]
                    errorApi=e 
                    ###########################
                allModel=dict()
                featureDict=dict()
                if sequenceAA=="" and sequenceNT=="":
                    sequenceAA="sequence not avalaible"
                    sequenceNT="sequence not avalaible"
                else:
                    ##########################feature engineering ####################################
                    data = {'geneSeq': sequenceNT, 'protSeq': sequenceAA}
                    args = {'geneSeq': sequenceNT, 'protSeq': sequenceAA}
                    try: 
                        r = requests.post(urlFeature, data=data, params=args, headers=header, verify=False)
                        feauturePredict=r.json()
                        CIARCSUFeat=feauturePredict['CIARCSUFeat']['DEG100011']
                        GCFeat=feauturePredict['GCFeat']['DEG100011']
                        GCFeat=GCFeat[0]
                        GeneLength=feauturePredict['GeneLength']['DEG100011']
                        KmerFeat=feauturePredict['KmerFeat']['DEG100011']
                        ProteinFeat=feauturePredict['ProteinFeat']['DEG100011']
                        lineFeatcsv=str(GeneLength)+","+str(GCFeat)+","+",".join([str(_) for _ in CIARCSUFeat])+","+",".join([str(_) for _ in KmerFeat])+","+",".join([str(_) for _ in ProteinFeat])+"\n"
                        lineFeattxt=str(GeneLength)+"  "+str(GCFeat)+" "+" ".join([str(_) for _ in CIARCSUFeat])+" "+" ".join([str(_) for _ in KmerFeat])+"    "+" ".join([str(_) for _ in ProteinFeat])+"\n"
                        #write to file
                        fileFeatcsv.write(str(lineFeatcsv))
                        fileFeattxt.write(str(lineFeattxt))

                    except Exception as e: 
                        CIARCSUFeat=""
                        GCFeat=""
                        GeneLength=""
                        KmerFeat=""
                        ProteinFeat=""
                        fileFeatcsv.write(str(e))
                        fileFeattxt.write(str(e))
                    
                    for  modelSelect, valueModel in modelDict.items():
                        modelResult=dict()
                        #sequenceNT=sequenceNT.upper()
                        #sequenceAA=sequenceAA.upper()
                        data = {'geneSeq': sequenceNT, 'protSeq': sequenceAA,'model':modelSelect}
                        args = {'geneSeq': sequenceNT, 'protSeq': sequenceAA,'model':modelSelect}
                        
                        ###############################predict ##########################################
                        try: 
                            r = requests.post(urlPredict, data=data, params=args, headers=header, verify=False)
                            resultPredict=r.json()
                            EssentialGeneProbability=resultPredict['EssentialGeneProbability']
                            NonEssentialGeneProbability=resultPredict['NonEssentialGeneProbability']
                            decision=resultPredict['decision']
                        except Exception as e: 
                            EssentialGeneProbability=""
                            NonEssentialGeneProbability=""
                            decision=""

                        linetxt=str(numberLine)+"   "+str(gene_locus)+"   "+str(proteinNameList[index])
                        linecsv=str(numberLine)+","+str(gene_locus)+","+str(proteinNameList[index])
                        linecsv=linecsv+","+str(valueModel)+","+str(EssentialGeneProbability)+","+str(NonEssentialGeneProbability)+","+str(decision)+","+str(sequenceNT)+","+str(sequenceAA)+"\n"
                        linetxt=linetxt+"   "+str(valueModel)+" "+str(EssentialGeneProbability)+"   "+str(NonEssentialGeneProbability)+"    "+str(decision)+"   "+str(sequenceNT)+" "+str(sequenceAA)+"\n"
                        filePredictcsv.write(str(linecsv))
                        filePredicttxt.write(str(linetxt))

                        modelResult={
                                'CIARCSUFeat':CIARCSUFeat,
                                'GCFeat':GCFeat,
                                'GeneLength':GeneLength,
                                'KmerFeat':KmerFeat,
                                'ProteinFeat':ProteinFeat,
                                'EssentialGeneProbability':EssentialGeneProbability,
                                'NonEssentialGeneProbability':NonEssentialGeneProbability,
                                'sequenceAA':sequenceAA,
                                'sequenceNT':sequenceNT,
                                'decision':decision
                                }
                        featureDict={
                                'CIARCSUFeat':CIARCSUFeat,
                                'GCFeat':GCFeat,
                                'GeneLength':GeneLength,
                                'KmerFeat':KmerFeat,
                                'ProteinFeat':ProteinFeat
                                }
                        #valueModel #allModel[modelSelect]=modelResult
                        allModel[valueModel]=modelResult
                        numberLine=numberLine+1
                        ##############################################################################
                item={
                    'officialName':gene_locus,
                    'organism_code':organism_code,
                    'protein_name':proteinNameList[index],
                    'sequenceAA':sequenceAA,
                    'sequenceNT':sequenceNT,
                    'geneLocus': gene_locus,
                    'model':allModel,
                    'modelDict':modelDict,
                    'feature': featureDict,
                    'errorSeq':errorSeq,
                    'errorApi':errorApi 
                    }

                resultDict[index]=item
                index=index+1
            #close file
            fileFeatcsv.close()
            fileFeattxt.close()
            filePredictcsv.close()
            filePredicttxt.close()

    return render(request, 'blog/ajaxSearchListProtein.html', locals())  
#vue du formulaire recherche organisms
def searchorganism(request):
    #load dataset 
    
    urlaws="http://3.133.132.243:9082/deeply/getKeggSpecie"

    url="http://195.24.221.74:5001/searchkeggorganis"
    species=dict()
    try: 
        r = requests.get(URL_API_SPECIES_KEGG)
        species=r.json()
        allSpecies=species 
    except Exception as err:
        try: 
            r = requests.get(urlaws)
            species=r.json()
        except Exception as error: 
            print("error: " + str(err))

    #init dictionary of organism
    
    #organismDict=dict()
    #for key, oneSpecies in species.items():
    #    specieId=oneSpecies['species_id']
    #    specieName=oneSpecies['official_name']
    #    organismDict[specieId]=specieName
    #text="bonjour"

    # Quoiqu'il arrive, on affiche la page du formulaire.
    return render(request, 'blog/searchorganism.html', locals())
#ajax proteinOrganism
def ajaxsearchorganism(request):

    if request.method == 'POST':
        action = request.POST.get('action', None)
        if action =='step1': 
            organism = request.POST.get('organism', None)
            step=action
            actionDict={"feature":"Feature Engeneering","prediction":"Prediction"}
            organism_code = request.POST.get('organism_code', None)
            url="http://3.133.132.243:9082/deeply/getProteinListEpath"
            data = {'organism_code': organism_code}
            args = {'organism_code': organism_code}
            header = {'Content-type': 'application/x-www-form-urlencoded; charset=utf-8'}
            listProtein=dict()
            try:
                r =requests.get(URL_API_PROTEIN_EPATH, params=args)
                listProtein=r.json()
            except :
                try: 
                    r = requests.post(url, data=data, params=args, headers=header, verify=False)
                    listProtein=r.json()
                except Exception as err:
                    print (err)
          

    return render(request, 'blog/ajaxsearchorganism.html', locals())
#predict by Text Mining
def predicttextmining(request):
    return render(request, 'blog/predicttextmining.html', locals())
#load you own model
def loadmodel(request):
    
    if request.method == 'POST':
        upload=True
        file=request.FILES['filename']
        #file = request.file['filename']
        if file and allowed_file(file.name):
            filename = secure_filename(file.name)
            content=file.read()
            fileSave = open(os.path.join(UPLOAD_FOLDER, filename), "w")
            fileSave.write(str(content))
            fileSave.close()
            #recupération des données du formulaire
            model_code=file.name
            modelname = request.POST.get('modelname', None)
            feauturemodel = request.POST.get('feauturemodel', None)
            modeldoi = request.POST.get('modeldoi', None)
            author = request.POST.get('author', None)
            authoraffiliation = request.POST.get('affiliationauthor', None)
            modeldescription = request.POST.get('modeldescription', None)
            
            data = {'model_code': model_code,'modelname':modelname,'feauturemodel':feauturemodel,'modeldoi':modeldoi,'author':author,'authoraffiliation':authoraffiliation,'modeldescription':modeldescription}
            args = {'model_code': model_code,'modelname':modelname,'feauturemodel':feauturemodel,'modeldoi':modeldoi,'author':author,'authoraffiliation':authoraffiliation,'modeldescription':modeldescription}
            header = {'Content-type': 'application/x-www-form-urlencoded; charset=utf-8'}
            url="http://3.133.132.243:9082/deeply/insertModel"
            message=""
            #request to get list of species share the protein name
            try: 
                r = requests.post(url, data=data, params=args, headers=header, verify=False)
                dumpJson=r.json()
                message=dumpJson['insert']
            except Exception as err: 
                message=err
            #file.save(os.path.join(UPLOAD_FOLDER, filename))
    else:
        #load a all model
        upload=False

    modelList=dict()
    try: 
        url="http://3.133.132.243:9082/deeply/getModels"
        r = requests.get(url)
        modelList=r.json()
    except Exception as err: 
        modelList=dict()


        

    return render(request, 'blog/loadmodel.html', locals())

#ajax load model
def ajaxloadmodel(request):

    if request.method == 'POST':

        return render(request, 'blog/ajaxloadmodel.html', locals())
        
#diplay API page
def displayapi(request):
    return render(request, 'blog/displayapi.html', locals())

#diplay colab page
def displaycolab(request):
    return render(request, 'blog/displaycolab.html', locals())

#diplay tutorial page
def displaytutorial(request):
    return render(request, 'blog/displaytutorial.html', locals())

def home(request):
    """ Exemple de page non valide au niveau HTML pour que l'exemple soit concis """
    return HttpResponse("""
        <h1>Bienvenue sur mon blog !</h1>
        <p>Les crêpes bretonnes ça tue des mouettes en plein vol !</p>
    """)


def view_article(request, id_article):
    # Si l'ID est supérieur à 100, nous considérons que l'article n'existe pas
    if id_article > 100:
        raise Http404

    return HttpResponse('<h1>Mon article ici</h1>')

def date_actuelle(request):
    date=datetime.now()
    return render(request, 'blog/date.html', locals())
    #return render(request, 'blog/date.html', {'date': datetime.now()})


def addition(request, nombre1, nombre2):    
    total = nombre1 + nombre2

    # Retourne nombre1, nombre2 et la somme des deux au tpl
    return render(request, 'blog/addition.html', locals())

#ajax predict by textmining
def ajaxpedicttextmining(request):

    if request.method == 'POST':
        action = request.POST.get('action', None)
        if action =='step1':
            elementid = request.POST.get('elementid', None)
            database = request.POST.get('database', None)
            actionDict={"feature":"Feature Engeneering","prediction":"Prediction"}
            urlSpeciesKegg="http://3.133.132.243:9082/deeply/searchKeggOrganism"
            header = {'Content-type': 'application/x-www-form-urlencoded; charset=utf-8'}
            step=action
            listProtein=dict()
            if database == 'pmc':
                # load dataset
                url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?pmcids="+elementid
                # r = requests.post(url, data=data, params=args, headers=header, verify=False)
                r = requests.get(url)
                pmc = r.json()
                allPUBMEDPMC = pmc
                # get pmc
            elif database == 'pubmed':
                # load dataset
                url = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?pmids="+elementid
                # r = requests.post(url, data=data, params=args, headers=header, verify=False)
                r = requests.get(url)
                pubmed = r.json()
                allPUBMEDPMC = pubmed
                # Get Pubmed
            # Parcours du JSON
            Title = allPUBMEDPMC["passages"][0]["text"]
            if database == 'pubmed':
                Doi = ""
            else:
                Doi = allPUBMEDPMC["passages"][0]["infons"]["article-id_doi"]
            Authors = allPUBMEDPMC["passages"][0]["infons"]["authors"]
            Year = allPUBMEDPMC["passages"][0]["infons"]["year"]
            if database == 'pubmed':
                Articleid_pmid = ""
            else:
                Articleid_pmid = allPUBMEDPMC["passages"][0]["infons"]["article-id_pmid"]
            if database == 'pubmed':
                Articleid_pmc = ""
            else:
                Articleid_pmc = allPUBMEDPMC["passages"][0]["infons"]["article-id_pmc"]
            index = 0
            index1 = 0
            tabresult = dict()
            speciestab = dict()

            for value in allPUBMEDPMC["passages"]:
                result = dict()
                result1 = dict()
                gene = ""
                species = ""
                taxaid = ""
                for value1 in value["annotations"]:
                    if value1["infons"]["type"] == "Gene":
                        gene = value1["text"]
                    if value1["infons"]["type"] == "Species":
                        species = value1["text"]
                        taxaid = value1["infons"]["identifier"]
                    if gene != "" and species != "":
                        result['gene'] = gene
                        result['species'] = species
                        result['taxaid'] = taxaid
                        result1['species'] = species
                        # test de l'existance de l'élément
                        verif = 0
                        if index > 0:
                            for test in tabresult:
                                if tabresult[test]['gene'] == gene and tabresult[test]['species'] == species:
                                    verif = 1
                        if verif == 0:
                            tabresult[index] = result
                            index = index + 1
                        
                        verif1 = 0
                        if index1 > 0:
                            for value1 in speciestab:
                                if speciestab[value1]['species'] == species:
                                    verif1 = 1
                        if verif1 == 0:
                            speciestab[index1] = result1
                            index1 = index1 + 1    
        elif action == 'step2':
            step2 = action
            # PREDICTION_VALUE = dict()
            # gene = request.POST.get('gene', None)
            speciesgenes = request.POST.get('speciesgenes', None)
            tabspeciesgenes = speciesgenes.split("/")
            genespeciesresult = dict()
            indexresult = 0
            for elt in tabspeciesgenes:
                smallresult = dict()
                smalltab = elt.split("_")
                smallresult = {"gene":smalltab[0],"species":smalltab[1],"taxaid":smalltab[2]}
                genespeciesresult[indexresult] = smallresult
                indexresult = indexresult + 1

            PREDICTION_VALUE = genespeciesresult

    return render(request, 'blog/ajaxpedicttextmining.html', locals())
