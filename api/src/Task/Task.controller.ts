import 'multer'

import {
  Body,
  Controller,
  Get,
  Inject,
  Post,
  UploadedFiles,
} from '@nestjs/common'

import { ClientRMQ } from '@nestjs/microservices'

import { FileStoreService } from './FileStore.service'
import { TaskService } from './Task.service'

import { UploadUnique } from '../common/UploadUnique'

export interface Task {
  id: string // TODO: make type-safe, share types with frontend
}

export interface UploadRequestBody {
  taskId: string
}

export interface PostTaskRequest {
  payload: {
    task: Task
  }
}

export interface GetTaskIdResponse {
  payload: {
    taskId: string // TODO: make type-safe, share types with frontend
  }
}

// TODO: return values should be different for different endpoints
export interface PostTaskResponse {
  payload: object
}

export interface UploadRequestFiles {
  DATES: Express.Multer.File[]
  FASTA: Express.Multer.File[]
  NWK: Express.Multer.File[]
}

@Controller()
export class TaskController {
  public constructor(
    private readonly taskService: TaskService,
    private readonly fileStoreService: FileStoreService,
    @Inject('TASK_QUEUE') private readonly taskQueue: ClientRMQ,
  ) {}

  @Get('/api/v1/taskId')
  public async getTaskId(): Promise<GetTaskIdResponse> {
    const taskId = this.taskService.generateTaskId()
    return { payload: { taskId } }
  }

  @Post('/api/v1/upload')
  @UploadUnique(['DATES', 'FASTA', 'NWK'])
  public async upload(
    @Body() { taskId }: UploadRequestBody,
    @UploadedFiles() files: UploadRequestFiles,
  ): Promise<PostTaskResponse> {
    const dates = files?.DATES?.[0]
    const fasta = files?.FASTA?.[0]
    const nwk = files?.NWK?.[0]

    const filesActual = [dates, fasta, nwk].filter(file => !!file)
    await this.fileStoreService.uploadInputFiles(taskId, filesActual)

    return { payload: { taskId } }
  }

  @Post('/api/v1/task')
  public async postTask(
    @Body() { payload: { task } }: PostTaskRequest,
  ): Promise<PostTaskResponse> {
    this.taskQueue.emit('tasks', 'hello')

    return { payload: { taskId: task.id } }
  }
}
