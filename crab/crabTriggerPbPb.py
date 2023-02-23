from WMCore.Configuration import Configuration
config = Configuration()

card='cardTriggerPbPb.input'
jobTag='triggerAnalysis_akFlowJets_eta1v6_baseCalo60_2023-02-22'
inputList='PbPb2018MiniAODTriggerEfficiency.txt'
outputFile=jobTag+'.root'
fileLocation='2'  # Locations: 0 = Purdue, 1 = CERN, 2 = Vanderbilt, 3 = Search with xrootd

config.section_("General")
config.General.requestName = jobTag
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+outputFile,'location='+fileLocation]
config.JobType.inputFiles = ['FrameworkJobReport.xml','jetTriggerAnalysis.tar.gz',card]
config.JobType.outputFiles = [outputFile]
config.JobType.maxJobRuntimeMin = 100
config.JobType.maxMemoryMB = 600

config.section_("Data")
config.Data.userInputFiles = open(inputList).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 100
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'PbPbTriggerHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_Vanderbilt']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

