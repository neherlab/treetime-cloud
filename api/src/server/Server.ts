import http from 'http'

import express from 'express'
import { useExpressServer } from 'routing-controllers'

import { getenv } from '../../lib/getenv'
import pkg from '../../package.json'
import controllers from '../controllers'
import { errorLogger, requestLogger } from './loggers'

export interface AppParams {
  useRequestLogger?: boolean
  useErrorLogger?: boolean
}

const appParamsDefault: Required<AppParams> = {
  useRequestLogger: true,
  useErrorLogger: true,
}

class Server {
  private readonly app: express.Application
  private readonly server: http.Server

  public constructor(params?: AppParams) {
    const { useRequestLogger, useErrorLogger } = {
      ...appParamsDefault,
      ...params,
    }

    this.app = express()
    this.app.set('etag', false)
    this.app.set('query parser', true)
    this.app.set('trust proxy', 'loopback')
    this.app.set('x-powered-by', false)

    this.server = http.createServer(this.app)

    if (useRequestLogger) {
      this.app.use(requestLogger())
    }

    useExpressServer(this.app, {
      cors: process.env.NODE_ENV === 'development',
      controllers,
    })

    if (useErrorLogger) {
      this.app.use(errorLogger())
    }
  }

  public get = () => {
    return this.app
  }

  public listen = () => {
    const port = getenv('API_PORT_INTERNAL')
    this.server.listen(port, () => {
      console.info(`${pkg.name} is listening on port ${port}`)
    })
  }

  public close = () => {
    this.server.close()
  }
}

export default Server
