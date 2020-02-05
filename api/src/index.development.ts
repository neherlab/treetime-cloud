import '../config/dotenv'

import { NestFactory } from '@nestjs/core'
import { NestExpressApplication } from '@nestjs/platform-express'

import { AppModule } from './app.module'
import { requestLogger } from './server/loggers'

import { getenv } from '../lib/getenv'
import pkg from '../package.json'

declare const module: NodeHotModule

async function bootstrap() {
  const app = await NestFactory.create<NestExpressApplication>(AppModule, {
    logger: ['error', 'warn'],
  })

  app.set('etag', false)
  app.set('query parser', true)
  app.set('trust proxy', 'loopback')
  app.set('x-powered-by', false)

  app.use(requestLogger())

  const port = getenv('API_PORT_INTERNAL')
  await app.listen(port)
  console.info(`${pkg.name} is listening on port ${port}`)

  // This will selectively hot-reload only the required parts of the application
  // on file changes in development mode
  if (module.hot) {
    module.hot.accept()
    module.hot.dispose(() => app.close())
  }
}

bootstrap()
