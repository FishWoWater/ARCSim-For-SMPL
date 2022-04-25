"""
This script runs a folder full of sims expecting/creating the following structure:
  folder/conf/*.json    ..  input
  folder/sims/          ..  output
"""
import platform
assert(platform.python_version().startswith("3")) # cluster safety
import os, sys, subprocess, argparse, time, datetime
import numpy as np
from multiprocessing import Pool
import getpass, smtplib, signal
from email.mime.text import MIMEText

email_user = ''
email_server = ''

# ARGS ########################################################################

ap = argparse.ArgumentParser()
ap.add_argument('folder', help='input folder with conf/ subfolder')

ap.add_argument("-b", "--build", default="1",
                help="...")
ap.add_argument("-r", "--run", default="1",
                help="...")
ap.add_argument("-o", "--output", default="sims",
                help="...")
ap.add_argument("-p", "--processes", default=1, type=int,
                help="...")
ap.add_argument("-t", "--threads", default=0, type=int,
                help="setting OMP_NUM_THREADS")
ap.add_argument("-d", "--delay", default=0, type=int,
                help="delay subprocess by n seconds to avoid initial remeshing clash")
ap.add_argument("-e", "--email", action='store_true')
ap.add_argument("-f", "--filter", default=[], nargs="*",
                help="filter to include")
ap.add_argument("-F", "--filterex", default=[], nargs="*",
                help="filter to exclude")
ap.add_argument("-q", "--queue", default=[], nargs="*",
                help="priority, e.g.: '-q basket 1 rib 2'")

args, unknownargs = ap.parse_known_args()
args = vars(args) # NOTE for future: instead of this can just use args.folder etc.
args['build'] = args['build'] != "0" and args['build'].lower() != "false"
args['run'] = args['run'] != "0" and args['run'].lower() != "false"


# HELPER FUNCS ################################################################

def try_timeout(func, T=60, default=None):
    signal.signal(signal.SIGALRM, lambda signum, frame: 1/0)
    signal.alarm(T)
    ret = default
    try:
        ret = func()
    except ZeroDivisionError:
        print("Timed out.")
        ret = default
    except:
        print("Aborted.")
        ret = default
    finally:
        signal.alarm(0)
    return ret

try_getmsg = lambda: try_timeout(lambda: input("Enter e-mail message: "))
try_getpass = lambda: try_timeout(lambda: getpass.getpass("Enter e-mail password: "))

def try_send_mail(header, content, server, user, pw):
    usermail = '%s@ist.ac.at' % user
    msg = MIMEText(content)
    msg['Subject'] = header
    msg['From'] = usermail
    msg['To'] = usermail

    try:
        server = smtplib.SMTP(server, 587)
        server.ehlo()
        server.starttls()
        server.login(user, pw)
        server.send_message(msg)
        print("Sent email.")
        server.close()
    except:
        print("Email error.")

def build(sourcedir, builddir, debug=False):
    cfg = 'Debug' if debug else 'Release'
    cmake_args = ['-DCMAKE_BUILD_TYPE=' + cfg]
    build_args = ['--config', cfg]
    build_args += ['--', '-j']

    os.makedirs(builddir, exist_ok=True)

    subprocess.check_call(['cmake', sourcedir] +
                          cmake_args, cwd=builddir)
    subprocess.check_call(['cmake', '--build', '.'] +
                          build_args, cwd=builddir)  # basically make

# BUILD #######################################################################

sourcedir = "."
builddir = os.path.join(sourcedir, "build-Release")

if args['build']:
    build(sourcedir=sourcedir, builddir=builddir, debug=False)

if not args['run']:
    exit()

env = {"OMP_NUM_THREADS": str(args['threads'])} if args['threads'] > 0 else None

if args['email']:
    email_pw = try_getpass()
    if email_pw is None or email_pw == "":
        args['email'] = False
    else:
        email_pw = email_pw[::-1]
        email_msg = try_getmsg()
        if email_msg is None:
            email_msg = ""

# RUN #########################################################################

workdir = sourcedir
executable = os.path.join(builddir, "bin", "arcsim_0.2.1")
# from workdir call: ./build-Release/bin/arcsim_0.2.1 $op $configfile ${3}

op = "simulateoffline"

# NOTE: no input sanitization, expecting -q expr1 priority1 expr2 priority2
confs = sorted(os.listdir(os.path.join(args['folder'], "conf"))) # default alphabetical
if len(args['queue']) > 0:
    priority_map = {args['queue'][2*i]: float(args['queue'][2*i+1]) for i in range(len(args['queue'])//2)}
    weights = [-sum([fil[1] for fil in priority_map.items() if fil[0] in conf]) for conf in confs]
    confs = np.array(confs)[np.argsort(weights, kind='mergesort')]

os.makedirs(os.path.join(args['folder'],args['output']), exist_ok=True)

confs_filtered = []
for i in range(len(confs)):
    # exclusive filter
    if any([fil in confs[i] for fil in args['filterex']]):
        continue
    # inclusive filter
    if not all([fil in confs[i] for fil in args['filter']]):
        continue
    confs_filtered.append(confs[i])
confs = confs_filtered

print("(Sorted/Filtered) Folders:\n" + "\n".join(confs))

t0 = datetime.datetime.now()

def tasks():
    delay = 0
    i = 0
    for conf in confs:
        if not conf.endswith(".json"):
            continue

        simname = os.path.splitext(os.path.basename(conf))[0]
        confpath = os.path.join(args['folder'], "conf", conf) 
        outputdir = os.path.join(args['folder'], args['output'], simname) 
        if os.path.isdir(outputdir):
            print("Skipping existing/in-progress", conf)
            continue

        yield (i,[executable, op, confpath, outputdir], delay)
        delay += args['delay']
        i += 1

def execute_task(args):
    i,task,delay = args
    outputdir = task[-1]
    if os.path.isdir(outputdir): # try skip again, if tasks preallocated then overlaps can happen..?
        print("Skipping existing/in-progress", task[-2])
        return
    if delay > 0:
        print("Delaying task %d for %02dm %02ds" % (i, delay//60, delay%60))
        time.sleep(delay)
    print("Starting task", i)
    subprocess.check_call(task, cwd=workdir, env=env)

try:
    n_processes = int(args['processes'])
    if n_processes > 1:
        pool = Pool(n_processes)
        iterator = tasks()
        pool_it = pool.imap_unordered(execute_task, iterator, chunksize=1)
        # wait and iterate pool results, keyboard abortable
        try:
            n_finished = 0
            for _ in pool_it:
                n_finished += 1
                print("Finished simulation #%d." % n_finished)
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate() # terminating 
        else:
            pool.close() # normal termination
    else:
        for task in tasks():
            execute_task(task)
    
except KeyboardInterrupt:
    print("PY: Aborting execution")

    if args['email']:
        print("Email message was:\n", email_msg)
        print("Elapsed:", datetime.datetime.now()-t0)

else:
    t1 = datetime.datetime.now()
    startstr   = "Started:  %d/%02d/%02d-%02d:%02d" % (t0.year, t0.month, t0.day, t0.hour, t0.minute)
    endstr     = "Finished: %d/%02d/%02d-%02d:%02d" % (t1.year, t1.month, t1.day, t1.hour, t1.minute)
    elapsedstr = "Elapsed:  %s" % (str(t1-t0))
    print(startstr, endstr, elapsedstr)

    if args['email']:
        try_send_mail("HYLC Finished", "The simulations have finished.\n" + email_msg + "\n\n%s\n%s\n%s" % (startstr,endstr,elapsedstr), email_server, email_user, email_pw[::-1])
