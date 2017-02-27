#!/usr/bin/env python
#786 

import os, sys, errno, argparse, subprocess, fnmatch, ConfigParser, time, stat, datetime

#from __future__ import print_function
#############################################################################################
# Class for colored texts and binary path
class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'

class pipeline:
	mechrna = os.path.dirname(os.path.realpath(__file__)) + "/mechrna.py"
	graphprot   = os.path.dirname(os.path.realpath(__file__)) + "/bin/graphprot_profile_prediction.pl"
	intarna  = os.path.dirname(os.path.realpath(__file__)) + "/bin/IntaRNA"
	intarna_pvalues  = "Rscript " + os.path.dirname(os.path.realpath(__file__)) + "/bin/intarna_pvals_gamma.R"
	inference  = "python " + os.path.dirname(os.path.realpath(__file__)) + "/bin/infer.py"
	workdir  = os.path.dirname(os.path.realpath(__file__))
	# example usage for help
	example  = "\tTo create a new project: specify (1) project name and (2) lncRNA sequence\n"
	example += "\t$ ./mechrna.py -p my_project -l my_rna.fa\n"	
	example += "\n\n\tTo resume a project, just type project folder and mechrna.py will automatically resume from the previous stages:\n"
	example += "\t$ ./mechrna.py -p my_project\n"
	example += "\t$ ./mechrna.py -p /home/this/is/my/folder/project\n\n"


#############################################################################################
# Default values will be set up later in check_proj_preq
def command_line_process():
	parser = argparse.ArgumentParser(
		description='MechRNA: LncRNA Mechanism Inference Tool',
		usage = pipeline.example,
		formatter_class=argparse.RawTextHelpFormatter
	)
	parser.add_argument('--project','-p',
		required=True,
		metavar='project',
		help='The name of the project. MechRNA creates the folder if it does not exist'
	)
	parser.add_argument('--lncRNA', '-l',
		required=True,
		metavar='lncRNA',
		help='Fasta file of lncRNA of interest.'
	)
	parser.add_argument('--rbps', '-b',
		metavar='rbps',
		help='List of RNA binding proteins to include in the analysis.'
	)
	parser.add_argument('--dest','-d',
		metavar='destination',
		help='Directory that will be used for analysis. (default: ./)'
	)
	parser.add_argument('--reference','-r',
		metavar='reference',
		help='"GRCh37.75" or "GRCh38.86"'
	)
	parser.add_argument('--topX','-x',
		type=int,
		metavar='topX',
		help="Top x percent of IntaRNA predictions to retain for further analysis. (default: 2)",
	)
	parser.add_argument('--cancer','-t',
		metavar='cancer',
		help="Cancer type of interest. Supported: prostate, meso.peri (default: prostate)",
	)
	parser.add_argument('--peak-cutoff','-P',
		type=float,
		metavar='peak_cutoff',
		help="P-value cutoff for protein binding peaks. (default: 0.01)",
	)
	parser.add_argument('--max-loop','-L',
		type=int,
		metavar='max_loop',
		help="Maximal distance of two paired bases in the local target RNA folding for computation of accessibilities (50nt < window). (default: 100)",
	)
	parser.add_argument('--window','-w',
		type=int,
		metavar='window',
		help="Size of the averaging window in the local target RNA folding for the computation of accessibilities. (default: 150)",
	)
	parser.add_argument('--suboptimals','-s',
		type=int,
		metavar='suboptimals',
		help="Number of suboptimals to compute. (default: 4)",
	)
	parser.add_argument('--max-len','-c',
		type=int,
		metavar='max_len',
		help="Maximum accessible transcript region. (default: 1000)",
	)
	parser.add_argument('--resume',
		nargs='?',
		const="intarna",
		help='Ignore existing progress and restart pipeline. Put intarna if you want to automatically resume from an previously killed task.',
	)
	parser.add_argument('--resume-force', '-f',
		action='store_true',
		help='Ignore existing files and restart the pipeline from stage specified in --resume.',
	)
	parser.add_argument('--num-worker',
		type=int,
		help='Number of independent prediction jobs which will be created. (default: 1)',
	)
	parser.add_argument('--range',
		help='Intervals of sequences in targets file to be analyzed.',
	)
	parser.add_argument('--worker-id',
		help='Specific worker ID that the user wants to run prediction in the last stage. (default:-1, will run on all workers)',
	)
	parser.add_argument('--mode',
		metavar='engine_mode',
		help='Type for running indepedent jobs: normal, sge, or pbs. (default: normal).',
		default='normal'
	)
	parser.add_argument('--job-num-cpus',
		help='Number of cpus per job in PBS file. Read documents before changing its value! (default: 1 cpu.)',
		default='1'
	)
	parser.add_argument('--job-max-time',
		help='Max job running time in PBS file. Read documents before changing its value! (default: 24 hour.)',
		default='24:00:00'
	)
	parser.add_argument('--job-max-memory',
		help='Max job memory in PBS file. Read docuemnts before changing its value! (default: 8 GB.)',
		default='8G'
	)

	return parser

#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def log( msg, output=True ):
	if output:
		print "{0}".format(msg),
		sys.stdout.flush()
	
#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def logln( msg, output=True ):
	if output:
		print "{0}".format(msg)
		sys.stdout.flush()

#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def logOK():
	logln(bcolors.OKGREEN+"OK"+bcolors.ENDC)

#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def logWAR():
	logln(bcolors.WARNING+"SKIPPING"+bcolors.ENDC)

#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def logFAIL():
	logln(bcolors.FAIL+"FAILED"+bcolors.ENDC)

#############################################################################################
### get the system shell to run the command and pipe stdout to control file when log==True
def shell(msg, run_command, command, control_file='', success_file='', success_msg='', shell_log=True):
	
	if shell_log:
		log("{0}...".format(msg))

	if (not run_command):
		logWAR()
		return 
	if not shell_log or control_file == '':
		p = subprocess.Popen(command, shell=True)
	else:
		f = open(control_file, 'w')
		f.write("CMD:" + command + "\n")
		f.flush() # so commands will show before actual progress
		p = subprocess.Popen(command, shell=True, stdout=f, stderr=f)

	ret = p.wait()
	if shell_log and control_file != '':
		f.close()
	if ret != 0:
		logFAIL()
		raise subprocess.CalledProcessError(ret, command)

	else:
		if shell_log and success_file != '':
			with open (success_file, 'w') as suc_file:
				suc_file.write(success_msg)
			logOK()

############################################################################################
##### check the status of a job. 0 :failed, 1: completed, 2: submitted either in the queue or running, 3: not submitted
def status_check(comp_file,success_message,pbs_err,pbs_log,PBS_TIME):
	st=3
	if not os.path.isfile(comp_file):
		st=3
	elif os.path.isfile(comp_file):
		comp_file_content =  open(comp_file,"r").read()
		if success_message in comp_file_content:
			st = 1
		elif 'submitted' in comp_file_content:
			st = 2
			if IsJobExpired(pbs_err,PBS_TIME):
				st=0
			elif os.path.isfile(pbs_err) and os.stat(pbs_err).st_size>0:
				err_mess = open(pbs_err,"r").read()
				if "job killed" in err_mess:
					if 'walltime' in err_mess:
						st=4
					elif 'nodes' in err_mess or 'ppn' in err_mess:
						st=5
					elif 'vmem' in err_mess:
						st=6
					elif 'Segmentation Fault' in err_mess:
						st=7
			elif os.path.isfile(pbs_log) and "FAILED" in open(pbs_log,"r").read():
				st=8
	return st
#############################################################################################
######### single 1 when it is a single job 0 when it is multiple jobs
def apply_status_check(st,single,message):
	if st == 0:
		logln("The job is already expired and there is no comp file")
		logFAIL()
		if single:
			exit(1)
	elif st == 1:
		logln(message+bcolors.WARNING+" SKIPPING"+bcolors.ENDC)
	elif st == 2:
		if single:
			logEXIT()
			exit(0)
	elif st==4:
		logln("Please check your pbs_err file, walltime might not be enough for this job.")
		if single:
			logEXIT()
			exit(0)
	elif st==5:
		logln("Please check your pbs_err file, nodes or ppn might not be enough for this job.")
		if single:
			logEXIT()
			exit(0)
	elif st==6:
		logln("Please check your pbs_err file, vmem might not be enough for this job.")
		if single:
			logEXIT()
			exit(0)
	elif st==7:
		logln("Please check your pbs_err file, out of memory.")
		if single:
			logEXIT()
			exit(0)
	elif st==8:
		logln("This is a cluster problem. Please check your pbs_err file. You can do nothing about this. Contact with cluster staff.\nOr you can try again:\n \
		In pbs_log delete the .o and .e of this job\n \
		In comp delete 00.calculateminmax.comp and 04.mappingsubmission.comp\n \
		In log delete 00.calculateminmax.log and 04.mappingsubmission.log\n \
		In mrsfastcomp delete .comp of this job\n \
		In mrsfastlog delete .log of this job\n \
		Decrease your submitted job numbers or threads in project.config file\n \
		Resubmit the job. It might still not work because it is a cluster error!.")
		if single:
			logEXIT()
			exit(0)
	return st
#############################################################################################
######### Return 1 if the error file exist longer than PBS running time. 0 otherwise: running
######### Absolute path of pbs.e file and a time string of format DD:MM:SS
def IsJobExpired( filename, PBS_TIME):	
	# The time since the pbs.error file is created
	if not os.path.isfile(filename):
		#print('.', end="")
		return 0
	cur=datetime.datetime.now()
	last=datetime.datetime.fromtimestamp(os.path.getmtime( filename ))
	tdelta = cur - last

	# The time of PBS time
	t_list = PBS_TIME.split(":")
	if 3 != len(t_list):
		print "Can not parse PBS TIME. Please supply in DD:MM:SS format"
		exit(1)
	
	pbs_hours = t_list[0]
	pbs_mins = t_list[1]
	pbs_secs = t_list[2]
	
	pbs_interval = datetime.timedelta(hours=int(pbs_hours), minutes=int(pbs_mins)+0, seconds=int(pbs_secs))

	flag = 0 # not exceed pbs time yet
	if (tdelta > pbs_interval):
		flag = 1
	return flag

#############################################################################################
### generate scripts to merge IntaRNA output when users run multiple workers
def generate_merge_script( config, num_worker ):
	appsdir  = os.path.dirname(os.path.realpath(__file__))
	worker_prefix   = "{0}/jobs/intarna".format( pipeline.workdir )
	output_prefix   = "{0}/{1}".format( pipeline.workdir, config.get("project", "name") )
	with open( pipeline.workdir + "/merge_output.sh", 'w' ) as f_script:
		f_script.write("#!/bin/bash\n")
		cmd = 'FILE={0}.results\nif [ -f $FILE ]; then rm -f ${{FILE}}; fi & touch ${{FILE}}\n'.format( output_prefix )
		f_script.write("#!/bin/bash\n{0}\n".format(cmd))
		for i in xrange( num_worker):
			f_script.write("part_info={0}_{1}; if [ -f ${{part_info}} ]; then cat ${{part_info}} >> ${{FILE}}; else printf \"Missing File %s\\n\" ${{part_info}};fi\n".format(worker_prefix, i))
	st = os.stat((pipeline.workdir + "/merge_output.sh"))
	os.chmod((pipeline.workdir + "/merge_output.sh"), st.st_mode | stat.S_IEXEC)

#############################################################################################
########## Clean stage file for resuming
def clean_state_worker( workdir, config):
	workdir = pipeline.workdir
	stage_dir = workdir + "/stage/"
	worker_stage_files = [f for f in os.listdir( stage_dir ) if (os.path.isfile( os.path.join( stage_dir , f)) and "02"==f[0:2])]
	for item in worker_stage_files:
		log(" + Removing " + item  + "...")
		os.remove( stage_dir  + item)
		logOK()
		
#############################################################################################
########## Clean stage file for resuming
def clean_state( mode_index, workdir, config ):
	flag_clean = 0 # 1 only when we delete files due to changes in project.config
	workdir = pipeline.workdir	
	valid_state = [ '01.graphprot', '02.worker', '04.pvalues', '05.inference', 'normal']
	for i in range( mode_index, len(valid_state)):
		if os.path.isfile(workdir + "/stage/" + valid_state[i] + ".finished" ):
		#try:
			if 0 == flag_clean:
				logln("Removing old files due to change in project.config")
				flag_clean = 1
			log(" + Removing " + valid_state[i] + ".finished...")
			os.remove( workdir + "/stage/" + valid_state[i] + ".finished" )
			logOK()
	# Clean all possible stages files of previous worker
	clean_state_worker( workdir, config)
	#except subprocess.CalledProcessError as e:
	#	print >>sys.stderr, "{0} failed with exit status {2} and message {1}".format(e.cmd, 'N/A', e.returncode)


#############################################################################################
###### Running commands for graphprot 
def graphprot(config ):
	msg           = "Predicting Protien Binding to lncRNA"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	input_file    = "{0}/{1}".format(workdir, config.get("project","lncRNA"))
	model         = "{0}/{1}".format(workdir, config.get("graphprot","rbps"))
	output_dir   = "{0}/{1}-profile".format(workdir, config.get("project","lncRNA"))
	graphprot_path= os.path.dirname(os.path.realpath(__file__)) + "/bin" 
	models_path   = config.get("graphprot","model-dir")
	control_file  = "{0}/log/01.graphprot.log".format(workdir);
	complete_file = "{0}/stage/01.graphprot.finished".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.graphprot + ' -fasta {0} -model {1} -out {2} -graphprot-path {3} -model-path {4}'.format(input_file, model, output_dir, graphprot_path, models_path)
	run_cmd       = not (os.path.isfile(complete_file) ) #TODO check if lncRNA exists
	shell( msg, run_cmd , cmd, control_file, complete_file, freeze_arg)

#############################################################################################
###### Running commands for intarna 
def intarna(config ):
	msg           = "Predicting RNA-RNA Interactions for transcripts {0}".format(config.get("intarna", "range"))
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	worker_id     = config.get("intarna", "worker-id")
	rng           = config.get("intarna", "range")
	target_file	  = "{0}/{1}".format(workdir, config.get("intarna","targets"))
	lncRNA_file   = "{0}/{1}".format(workdir, config.get("project","lncRNA"))
	output_prefix   = "{0}/jobs/intarna_{1}".format(workdir, worker_id )
	control_file  = "{0}/log/03.intarna.{1}.log".format(workdir, worker_id);
	complete_file = "{0}/stage/03.intarna.{1}.finished".format(workdir, worker_id);
	success_message	 = "completed"
	#freeze_arg    = "{0}\t{1}\t{2}\t{3}\t{4}".format(worker_id, rng, config.get("intarna","max-loop"), config.get("intarna","window"), config.get("intarna","suboptimals"))
	cmd           = pipeline.intarna +' -P -L {0} -w {1} -s {2} -c {3} -t {4} -m {5} -r {6} > {7}'.format(config.get("intarna","max-loop"), config.get("intarna","window"), config.get("intarna","suboptimals"), config.get("intarna","max-len"), target_file, lncRNA_file, rng, output_prefix)
	run_cmd       = not ( os.path.isfile(complete_file) and success_message in open(complete_file).read()) 
	shell( msg, run_cmd, cmd, control_file, complete_file, success_message)

#############################################################################################
###### Running commands for intarna pvalues 
def intarna_pvalues(config):
	msg           = "Computing P-values"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	input_file    = "{0}/{1}.results".format(workdir, config.get("project", "name"))
	control_file  = "{0}/log/04.pvalues.log".format(workdir);
	complete_file = "{0}/stage/04.pvalues.finished".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.intarna_pvalues + ' {0} {1}'.format(config.get("project", "topX"), input_file)
	run_cmd       = not ( os.path.isfile(complete_file) and freeze_arg in open(complete_file).read()) 

	if ( run_cmd ):
		clean_state( 2, workdir, config)
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
	
#############################################################################################
###### Running commands for inference
def inference(config):
	msg           = "Inferring Mechanisms"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	input_file    = "{0}/{1}.results.pvalues".format(workdir, config.get("project","name"))
	models        = "{0}/{1}".format(workdir, config.get("graphprot","rbps"))
	output_file   = "{0}/{1}.results.pvalues.mechs".format(workdir, config.get("project","name"))
	control_file  = "{0}/log/05.inference.log".format(workdir);
	complete_file = "{0}/stage/05.inference.finished".format(workdir);
	freeze_arg    = ""
	cmd           = pipeline.inference + ' -i {0} -r {1} -m {2} -c {3} -p {4} > {5}'.format(input_file, config.get("project", "reference"), models, config.get("project", "cancer"), config.get("project", "peak-cutoff"), output_file ) 
	run_cmd       = not ( os.path.isfile(complete_file) and freeze_arg in open(complete_file).read()) 

	if ( run_cmd ):
		clean_state( 3, workdir, config )
	shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)

#############################################################################################
########## Generate commands for each worker
def assign_worker(config):
	msg           = "Worker"
	project_name  = config.get("project", "name")
	workdir		  = pipeline.workdir
	target_file	  = "{0}/{1}".format(workdir, config.get("intarna","targets"))
	lncRNA_file   = "{0}/{1}".format(workdir, config.get("project","lncRNA"))
	output_file   = "{0}/jobs.sh".format(workdir)
	control_file  = "{0}/log/02.worker.log".format(workdir);
	complete_worker_file = "{0}/stage/02.worker.finished".format(workdir);
	nj            = 0;

	with open(target_file) as fasta:
		for line in fasta:
			if (line[0] == ">"):
				nj=nj+1

	maxjobs = int( config.get("project", "num-worker") )
	if (maxjobs < 1):# sanity checking
		maxjobs = 1 
	
	freeze_arg    = "{0}\t{1}\t{2}\t{3}\t{4}".format(nj, maxjobs, config.get("intarna","max-loop"), config.get("intarna","window"), config.get("intarna","suboptimals"))

	run_cmd       = not ( os.path.isfile(complete_worker_file) and freeze_arg in open(complete_worker_file).read()) 
	if (run_cmd):
		log("Assign jobs to independent workers...")
		with open (control_file, 'a') as ctr_file:
			ctr_file.write("Find {0} sequences\n".format(nj))
		#maxjobs = int( config.get("project", "num-worker") )
		#if (maxjobs < 1):# sanity checking
		#	maxjobs = 1 
	
		if os.path.isfile(control_file ):
			os.remove( control_file)

		#with open(output_file, 'w') as fj:
		#	fj.write("#!/bin/bash\n")
		
		for i in xrange(maxjobs):

			complete_file = "{0}/stage/03.intarna.{1}.finished".format(workdir, i);
			if( not os.path.isfile(complete_file)):
				# worker i
				st = (nj / maxjobs) * i
				if st >= nj:
					break
				ed = min(nj, ((nj / maxjobs) * (i + 1)-1))
				rng = '{0}-{1}'.format(st, ed)
				cmd = pipeline.mechrna + " -p {0} -l {1} --range {2} --worker-id {3} --num-worker 1 --resume intarna \n".format( pipeline.workdir, lncRNA_file, rng, i )
				if "sge" == config.get("intarna", "engine-mode"):
					filename = '{0}/pbs/intarna_{1}.sh'.format(workdir, i)
					with open(filename, 'w') as fsge:
						pbs_log=pipeline.workdir+'/log/intarna_{0}.o'.format(i)
						pbs_err= pipeline.workdir+'/log/intarna_{0}.e'.format(i)
						worker_cmd = "qsub -cwd -V -b y -N {0}_{1} -pe ncpus {2} -l h_vmem={3} -l h_rt={4} -l h_stack=8M -o {5} -e {6} ".format( project_name, i, config.get("intarna", "job-cpus"), config.get("intarna", "job-memory"), config.get("intarna","job-time"), pbs_log, pbs_err) +  cmd;
						fsge.write("#!/bin/bash\n")
						fsge.write(worker_cmd)
						st = os.stat(filename)
						os.chmod(filename, st.st_mode | stat.S_IEXEC)
				elif "pbs" == config.get("intarna", "engine-mode"):
					with open('{0}/pbs/intarna_{1}.pbs'.format(workdir, i), 'w') as fpbs:
						fpbs.write("#!/bin/bash\n")
						fpbs.write("#PBS -l nodes=1:ppn={0},vmem={1},walltime={2}\n".format(config.get("intarna", "job-cpus"), config.get("intarna", "job-memory"), config.get("intarna", "job-time")) ) 
						fpbs.write("#PBS -N {0}_{1}\n".format(project_name, i))
						fpbs.write("cd $PBS_O_WORKDIR\n")
						fpbs.write("{0}".format(cmd))
					worker_cmd = "qsub {0}/pbs/{1}.pbs\n".format(workdir,i)
				else: # local machine
					worker_cmd = cmd + "\n"
	
				#with open('{0}'.format(output_file), 'a') as fj:
				#	fj.write(worker_cmd)

				with open (control_file, 'a') as ctr_file:
					ctr_file.write("Finish worker {0}\n".format(i))
				#log("Job {0}...".format(i))
				#logOK()
			else:
				log("Job {0}...".format(i))
				logWAR()

		#st = os.stat(output_file)
		#os.chmod(output_file, st.st_mode | stat.S_IEXEC)

		generate_merge_script( config, maxjobs  )
		with open (complete_worker_file, 'w') as suc_file:
			suc_file.write("{0}".format( freeze_arg ))
			logOK()
		# to re-run jobs
		clean_state( 10, workdir, config )

###############################################################################################			
########### submitting IntaRNA jobs
def submit_intarna_jobs(config):
	maxjobs = int( config.get("project", "num-worker") )
	workdir	= pipeline.workdir
	success_message	  = "completed"
	submitted = False
	failed=0
	running=0
	engine_mode = config.get("intarna", "engine-mode")
	for i in xrange(maxjobs):
		if("pbs" ==  engine_mode):
			worker_cmd = "{0}/pbs/intarna_{1}.pbs".format(workdir,i)
		elif("sge" == engine_mode):
			worker_cmd = "{0}/pbs/intarna_{1}.sh".format(workdir,i)
		pbs_log=pipeline.workdir+'/log/intarna_{0}.o'.format(i)
		pbs_err= pipeline.workdir+'/log/intarna_{0}.e'.format(i)
		complete_file = "{0}/stage/03.intarna.{1}.finished".format(workdir, i);
		st=status_check(complete_file,success_message,pbs_err,pbs_log,config.get("intarna", "job-time"))
		st=apply_status_check(st,False,"Submitting PBS job for worker {0}".format(i))
		if st == 0:
			failed+=1
		elif st == 2:
			running+=1
		elif st == 3:
			if("pbs" ==  engine_mode):
				shell(" + Submitting intarna_{0}.pbs".format(i),True, "qsub "+worker_cmd+" -o "+pbs_log+' -e '+pbs_err,'',complete_file,'submitted')
			elif("sge" == engine_mode):
				shell(" + Submitting intarna_{0}.sh".format(i),True, worker_cmd,'',complete_file,'submitted')
			submitted=True

	return running,failed

#############################################################################################
###### Running commands for each mode 
def run_command(config, force=False):

	#graphprot(config) #TODO run this if transcript is novel, add check
	assign_worker(config)
	success_message	  = "completed"

	engine_mode = config.get("intarna", "engine-mode")
	if ( ( "pbs" !=  engine_mode) and ( ( "sge" ) != engine_mode) ):
		num_parallel = 1
		msg = "Running prediction with {0} parallel tasks".format( num_parallel )
		script_file =  pipeline.workdir+ "/jobs.sh"
		cmd = "cat {0} | xargs -I CMD --max-procs={1} bash -c CMD ".format(script_file, num_parallel )
		freeze_arg=""
		control_file  = "{0}/log/normal.log".format( pipeline.workdir )
		complete_file = "{0}/stage/normal.finished".format( pipeline.workdir)
		run_cmd       = not ( os.path.isfile(complete_file) and freeze_arg in open(complete_file).read()) 
		shell( msg, run_cmd, cmd, control_file, complete_file, freeze_arg)
	else:
		#msg = "Submitting jobs to cluster\n")
		#cmd = "{0}/jobs.sh".format(pipeline.workdir)
		#shell(msg, True, cmd)
		submit_intarna_jobs(config)
		logOK()

	jobs_finished = False;
	maxjobs = int( config.get("project", "num-worker") )
	if (maxjobs < 1):# sanity checking
		maxjobs = 1 

	while(not jobs_finished):
		time.sleep(60)
		jobs_finished = True
		print ".",
		for i in xrange(maxjobs):
			complete_file = "{0}/stage/03.intarna.{1}.finished".format(pipeline.workdir, i);
			pbs_log=pipeline.workdir+'/log/intarna_{0}.o'.format(i)
			pbs_err= pipeline.workdir+'/log/intarna_{0}.e'.format(i)

			st=status_check(complete_file,success_message,pbs_err,pbs_log,config.get("intarna", "job-time"))
			st=apply_status_check(st,False,"")

			if not st == 1:
				jobs_finished = False

	msg = "Merging {0} output files".format(config.get("project", "num-worker"))
	cmd = "{0}/merge_output.sh".format(pipeline.workdir)
	shell( msg, True, cmd)
	logOK()

	intarna_pvalues(config)
	inference(config)

#############################################################################################
def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as e:
		if e.errno == errno.EEXIST and os.path.isdir(path):
			print "[ERROR] The project folder exists. Please run in resume mode or delete the project folder to re-run from scartch"
			exit(1);

#############################################################################################
######### link the absolute path of src in dest with identical filename
def symlink(src, dest):
	if src == '':
		return
	if not os.path.isfile(src):
		print "[ERROR] Input file {0} does not exist".format(src)
		exit(1)
	basename = os.path.basename(src)
	dest=os.path.abspath(dest)+'/'+basename
	#abs_src= os.path.abspath(src) + "/" + basename
	#os.symlink(src, dest)
	os.symlink(os.path.abspath(src), dest)

#############################################################################################
######### link the absolute path of src in dest with new filename
def symlink_name(src, dest, filename ):
	if src == '':
		return
	if not os.path.isfile(src):
		print "[ERROR] Input file {0} does not exist".format(src)
		exit(1)
	basename = os.path.basename(filename)
	dest=os.path.abspath(dest)+'/'+basename
	#abs_src= os.path.abspath(src) + "/" + basename
	#os.symlink(src, dest)
	os.symlink(os.path.abspath(src), dest)
#############################################################################################
##########	Make sure the project is successfully built
def is_exec(f):
	return os.path.isfile(f) and os.access(f, os.X_OK)

def check_binary_preq():
	execs = ['IntaRNA', 'Rscript']#, 'perl']
	log( "Checking binary pre-requisites... ")
	for exe in execs:
		local_exe = os.path.dirname(os.path.realpath(__file__)) + "/bin/" + exe
		installed = False
		if not is_exec(local_exe):
			for path in os.environ["PATH"].split(os.pathsep):
				path = path.strip('"')
				installed_exe = os.path.join(path, exe)
				if(is_exec(installed_exe)):
					installed = True
		else:
			installed = True

		if(not installed):
			if(exe == 'IntaRNA'):
				print "[ERROR] File {0} cannot be executed. Please use 'make' to build the required binaries.".format(exe)
				logFAIL()
				logln ("File {0} cannot be executed. Please use 'make' to build the required binaries.".format(exe) )
			else:
				print "[ERROR] File {0} is not installed in system PATH or in {1}".format(exe, os.path.dirname(os.path.realpath(__file__)) + "/bin/")
				logFAIL()
				logln ("File {0} is not installed in system PATH or in {1}".format(exe, os.path.dirname(os.path.realpath(__file__)) + "/bin/") )
			exit(1)
#TODO check for R packages: evd, 

	logOK()

#############################################################################################
def resume_state_help():
	print "\nMechRNA supports the following resume states:"
	print "\tgraphprot: computes protein binding sites on the lncRNA"
	print "\tintarna: predicts RNA-RNA interactions"
	print "\tpvalues: computes intarna p-values"
	print "\tinference: determines potential mechanism using correlation and relative locations of interactions"
	print "\tnum-worker: generate jobs for parallel processing"
	print "\nNOTE\tIf you want to automatically resume a killed job, just type --resume"

#############################################################################################
# Checking if necessary files are provided for a NEW project.
def check_input_preq( config ):
	workdir = pipeline.workdir

	if (None == config.get("project", "lncRNA") ) or (not os.path.isfile( config.get("project", "lncRNA") ) ):
		logFAIL()
		logln("The lncRNA sequence file, {0}, does not exist. Please provide valid name or path.".format( config.get("project", "lncRNA") ))
		exit(1)

#############################################################################################
########## Initialze mrsfast parameters for before creating project folder
def initialize_config_graphprot( config, args):
	config.add_section("graphprot")
	config.set("graphprot", "rbps", str( args.rbps ) if args.rbps != None else  "{0}/data/models.list".format(os.path.dirname(os.path.realpath(__file__))))
	config.set("graphprot", "model-dir", "{0}/data/models".format(os.path.dirname(os.path.realpath(__file__))))
	return config #TODO structue vs sequnence model

#############################################################################################
########## Initialze intarna parameters for before creating project folder
def initialize_config_intarna( config, args):
	workdir = pipeline.workdir
	config.add_section("intarna")
	config.set("intarna", "max-loop", str( args.max_loop ) if args.max_loop != None else  "150")
	config.set("intarna", "window", str( args.window ) if args.window != None else "200")
	config.set("intarna", "suboptimals", str( args.suboptimals ) if args.suboptimals !=None else "4")
	config.set("intarna", "max-len", str( args.max_len ) if args.max_len !=None else "1000")
	config.set("intarna", "targets", "{0}/data/refs/ensembl.{1}.ncrna.cdna.40".format(os.path.dirname(os.path.realpath(__file__)), config.get("project", "reference")))
#	config.set("intarna", "targets", "{0}/data/test.targets.fa".format(os.path.dirname(os.path.realpath(__file__))))
	config.set("intarna","engine-mode",args.mode if args.mode != None else "normal")
	config.set("intarna","range", str( args.range ) if args.range !=None else "-1" )
	config.set("intarna","worker-id", str(args.worker_id) if args.worker_id != None else "-1")
	config.set("intarna","job-cpus", args.job_num_cpus if args.job_num_cpus != None else "1")
	config.set("intarna","job-time", args.job_max_time if args.job_max_time != None else "24:00:00")
	config.set("intarna","job-memory",args.job_max_memory if args.job_max_memory != None else "8G")
	return config

#############################################################################################
def check_project_preq():
	args = command_line_process().parse_args()
	config = ConfigParser.ConfigParser()
	
	project_name = os.path.basename(os.path.normpath(args.project))
	pipeline.workdir = os.path.abspath(args.project)
	workdir = pipeline.workdir

	print "============================================="
	print "Project Name      : "+bcolors.OKGREEN+project_name+bcolors.ENDC
	print "Working Directory : "+bcolors.OKGREEN+workdir+bcolors.ENDC
	print "============================================="
	
	# Check if users want to resume the project
	if ( os.path.isdir( workdir)):
		log ("Checking the project pre-requisites... ")
		if ( None == args.resume):
			logFAIL()
			logln("MechRNA can not overwrite an existing project. Please add --resume or change project name.")
			exit(1)
			
		if not os.path.isfile( workdir + "/project.config"):
			logFAIL()
			logln("NO config settings found. Please remove and re-create the project.")
			exit(1)
		logOK()
		
		log ("Loading the config file... ")
		config.read(workdir + '/project.config');
		# update range  and worker id for assemble stage in SGE and PBS
		config.set("intarna","range", str( args.range ) if args.range !=None else "-1" )
		config.set("intarna","worker-id", str( args.worker_id ) if args.worker_id !=None else "-1" )
		logOK()
	
	# Start a new project
	else:
		log("Creating a new project folder...")

		if ((args.reference != None) and (args.reference != "GRCh37.75") and (args.reference != "GRCh38.86")):
			logFAIL()
			logln("Invalid reference specified. Please use \"GRCh37.75\" or \"GRCh38.86\"")
			exit(1)

		# set up main project parameter
		config.add_section("project")
		config.set("project", "name", project_name)
		config.set("project", "lncRNA", args.lncRNA)
		config.set("project", "reference", str(args.reference) if args.reference != None else "GRCh37.75") 
		config.set("project", "num-worker", str(args.num_worker) if args.num_worker != None else "1" )
		config.set("project", "topX", args.topx if args.topx != None else "2" )
		config.set("project", "cancer", str(args.cancer) if args.cancer != None else "prostate" )
		config.set("project", "peak-cutoff", args.peak_cutoff if args.peak_cutoff != None else "0.01" )

		# Parameters for other parts in the pipeline
		initialize_config_graphprot(config, args)
		initialize_config_intarna(config, args)
		
		#validating required files according to mode
		check_input_preq(config)
		
		# creating project folder
		mkdir_p(workdir)
		mkdir_p(workdir +'/jobs');
		mkdir_p(workdir +'/log');
		mkdir_p(workdir +'/pbs');
		mkdir_p(workdir +'/stage');

		symlink(config.get("graphprot","rbps"),  workdir)
		config.set("graphprot", "rbps", os.path.basename(config.get("graphprot", "rbps")))
		
		symlink(config.get("intarna", "targets"),  workdir)
		config.set("intarna", "targets", os.path.basename(config.get("intarna", "targets")))

		symlink(config.get("project", "lncRNA"),  workdir)
		config.set("project", "lncRNA", os.path.basename(config.get("project", "lncRNA")))

		# creating config in folder
		with open ( workdir +"/project.config", "w") as configFile:
			config.write(configFile)
		logOK()

	return config
#############################################################################################
def main():
	config = check_project_preq()
	check_binary_preq()

	resume_state="intarna"

	try:
		if config.get("intarna", "range") == '-1':
			run_command(config)
		elif resume_state == "intarna":
			intarna(config)
		else:
			raise Exception('Invalid mode selected: ' + mode)
	except subprocess.CalledProcessError as e:
		print >>sys.stderr, "{0} failed with exit status {2} and message {1}".format(e.cmd, 'N/A', e.returncode)

#############################################################################################
if __name__ == "__main__":
    sys.exit(main())

