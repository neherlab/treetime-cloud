### Webserver to run TreeTime analysis

#### Install

```shell
npm install
webpack
```
The latter build the java script file from the react code.

The python requirements are
  * TreeTime
  * Flask
  * werkzeug
  * threading

#### Run

```shell
python treetime_server.py
```
will fire up a local webserver on (default) port 4100.