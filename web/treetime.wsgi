#!/usr/bin/python
import sys
import logging

root_dir =  '/var/www/treetime_web'
treetime_dir =  '/var/www/treetime'


if not treetime_dir in sys.path:
	sys.path.insert(0, treetime_dir)

if not root_dir in sys.path:
	sys.path.insert(0,root_dir)

logging.basicConfig(stream=sys.stderr)

from treetime_server import app
application = app

