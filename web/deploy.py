"""
Web server deploy script
"""

import os
import subprocess


host = 'treetime_test'
host_location = '/var/www/treetime_web'
treetime_admin = 'root'
treetime_user = 'ttime'

if __name__ == '__main__':
    # first, copy all necessary files to the remote host
    os.system('scp -r ./examples/ ./static/ ./templates/ ./treetime.wsgi ./treetime_server.py {}@{}:{}/'.format(treetime_admin, host, host_location))

    # check if the sessions folder is in place
    if not 'sessions' in subprocess.check_output('ssh -t {}@{}  "ls {}"'.format(treetime_admin, host, host_location), shell=True):
        os.system('ssh -t {}@{}  "mkdir {}/sessions"'.format(treetime_admin, host, host_location))
        os.system('ssh -t {}@{}  "chown -R {}:{} {}/sessions"'.format(treetime_admin, host, treetime_user, treetime_user, host_location))

    #change permissions of the wsgi script:
    os.system('ssh -t {}@{}  "chown {}:{} {}/treetime.wsgi"'.format(treetime_admin, host, treetime_user, treetime_user, host_location))






