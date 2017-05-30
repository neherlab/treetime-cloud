"""
Web server deploy script
"""

import os

command = 'scp'
host = 'treetime_test'
host_location = '/var/www/treetime_web'


if __name__ == '__main__':

    #os.system('{} -r {}:{}/examples/ ./examples/'.format(command, host, host_location))
    #os.system('{} -r {}:{}/static/ ./static/'.format(command, host, host_location))
    #os.system('{} -r {}:{}/templates/ ./templates/'.format(command, host, host_location))
    os.system('{} ./treetime.wsgi {}:{}'.format(command, host, host_location))
    os.system('{} ./treetime_server.py {}:{}'.format(command, host, host_location))


