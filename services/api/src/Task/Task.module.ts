import { Inject, Module } from '@nestjs/common'
import { ClientRMQ, ClientsModule, Transport } from '@nestjs/microservices'

import { getenv } from '../../lib/getenv'

import { TaskController } from './Task.controller'
import { FileStoreService } from './FileStore.service'
import { TaskIdService } from './TaskId.service'

const TASK_QUEUE_HOST = getenv('TASK_QUEUE_HOST')
const TASK_QUEUE_PORT = getenv('TASK_QUEUE_PORT')
const TASK_QUEUE_ADDRESS = `${TASK_QUEUE_HOST}:${TASK_QUEUE_PORT}`

const RmqClientModule = ClientsModule.register([
  {
    name: 'TASK_QUEUE',
    transport: Transport.RMQ,
    options: {
      urls: [TASK_QUEUE_ADDRESS],
      queue: 'tasks',
      noAck: false,
      queueOptions: { durable: false },
    },
  },
])

@Module({
  imports: [RmqClientModule],
  controllers: [TaskController],
  providers: [TaskIdService, FileStoreService],
})
export class TaskModule {
  constructor(@Inject('TASK_QUEUE') private readonly taskQueue: ClientRMQ) {}

  private async onApplicationBootstrap() {
    await this.taskQueue.connect()
  }
}
