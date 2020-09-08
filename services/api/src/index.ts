import '../config/dotenv'

import { NestFactory } from '@nestjs/core'
import { NestExpressApplication } from '@nestjs/platform-express'
import { MicroserviceOptions, Transport } from '@nestjs/microservices'

import { AppModule } from './App.module'
import { requestLogger } from './common/logger.middleware'

import { getenv } from '../lib/getenv'
import pkg from '../package.json'

declare const module: NodeHotModule

const TASK_QUEUE_HOST = getenv('TASK_QUEUE_HOST')
const TASK_QUEUE_PORT = getenv('TASK_QUEUE_PORT')
const TASK_QUEUE_ADDRESS = `${TASK_QUEUE_HOST}:${TASK_QUEUE_PORT}`

console.log({ TASK_QUEUE_ADDRESS })

async function bootstrap() {
  const httpServer = await NestFactory.create<NestExpressApplication>(AppModule, {
    logger: ['error', 'warn'],
  })

  const rmqConsumer = await NestFactory.createMicroservice<MicroserviceOptions>(AppModule, {
    transport: Transport.RMQ,
    options: {
      urls: [TASK_QUEUE_ADDRESS],
      queue: 'taskResults',
      noAck: false,
      queueOptions: { durable: false },
    },
  })

  httpServer.set('etag', false)
  httpServer.set('query parser', true)
  httpServer.set('trust proxy', 'loopback')
  httpServer.set('x-powered-by', false)

  httpServer.use(requestLogger())

  const port = getenv('API_PORT_INTERNAL')

  await Promise.all([rmqConsumer.listenAsync(), httpServer.listenAsync(port)])
  console.info(`${pkg.name}: HTTPS server is listening on port ${port}`)
  console.info(`${pkg.name}: RabbitMQ consumer is listening`)

  // This will selectively hot-reload only the required parts of the application
  // on file changes in development mode
  if (module.hot) {
    module.hot.accept()
    module.hot.dispose(async () => {
      await httpServer.close()
      await rmqConsumer.close()
    })
  }
}

bootstrap().catch(console.error)
