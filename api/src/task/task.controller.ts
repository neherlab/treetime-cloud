import {
  Body,
  Controller,
  Get,
  Inject,
  Post,
  UploadedFiles,
  UseInterceptors,
} from '@nestjs/common'

import { FileFieldsInterceptor } from '@nestjs/platform-express'

import { ClientRMQ } from '@nestjs/microservices'

import { TaskService } from './task.service'

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
  DATES: File[]
  FASTA: File[]
  NWK: File[]
}

@Controller()
export class TaskController {
  // HACK: should become a service, with client isolation
  private files: Map<string, File> = new Map<string, File>()

  public constructor(
    private readonly taskService: TaskService,
    @Inject('TASK_QUEUE') private readonly taskQueue: ClientRMQ,
  ) {}

  @Get('/api/v1/taskId')
  public async getTaskId(): Promise<GetTaskIdResponse> {
    const taskId = this.taskService.generateTaskId()
    return { payload: { taskId } }
  }

  @Post('/api/v1/upload')
  @UseInterceptors(
    FileFieldsInterceptor([
      { name: 'DATES', maxCount: 1 },
      { name: 'FASTA', maxCount: 1 },
      { name: 'NWK', maxCount: 1 },
    ]),
  )
  public async upload(
    @Body() { taskId }: UploadRequestBody,
    @UploadedFiles() files: UploadRequestFiles,
  ): Promise<PostTaskResponse> {
    const dates = files?.DATES?.[0]
    const fasta = files?.FASTA?.[0]
    const nwk = files?.NWK?.[0]

    if (dates) {
      this.files.set('DATES', dates)
    }

    if (fasta) {
      this.files.set('FASTA', fasta)
    }

    if (nwk) {
      this.files.set('NWK', nwk)
    }

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
