import { Inject, Module } from '@nestjs/common'
import { ClientRMQ, ClientsModule, Transport } from '@nestjs/microservices'

import { TaskController } from './Task.controller'

import { FileStoreService } from './FileStore.service'
import { TaskIdService } from './TaskId.service'

const RmqClientModule = ClientsModule.register([
  {
    name: 'TASK_QUEUE',
    transport: Transport.RMQ,
    options: {
      urls: [`amqp://treetime-dev-taskqueue:5672`],
      queue: 'tasks',
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
