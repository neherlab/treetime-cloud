import 'multer'

import {
  BadRequestException,
  Body,
  Controller,
  Get,
  Inject,
  NotFoundException,
  Post,
  UploadedFiles,
} from '@nestjs/common'

import { ClientRMQ } from '@nestjs/microservices'

import serialize from 'serialize-javascript'

import { UploadUnique } from '../common/UploadUnique'

import { FileStoreService, UploadedFileData } from './FileStore.service'
import { TaskIdService } from './TaskId.service'

export interface Task {
  taskId: string // TODO: make type-safe, share types with frontend
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
  payload: Record<string, unknown>
}

@Controller()
export class TaskController {
  public constructor(
    private readonly taskIdService: TaskIdService,
    private readonly fileStoreService: FileStoreService,
    @Inject('TASK_QUEUE') private readonly taskQueue: ClientRMQ,
  ) {}

  @Get('/api/v1/taskId')
  public async getTaskId(): Promise<GetTaskIdResponse> {
    const taskId = this.taskIdService.generateTaskId()
    return { payload: { taskId } }
  }

  @Post('/api/v1/upload')
  @UploadUnique(['DATES', 'FASTA', 'NWK'])
  public async upload(
    @Body() { taskId }: UploadRequestBody,
    @UploadedFiles() filedata: UploadedFileData,
  ): Promise<PostTaskResponse> {
    await this.fileStoreService.uploadInputFiles(taskId, filedata)
    return { payload: { taskId } }
  }

  @Post('/api/v1/task')
  public async postTask(@Body() body: PostTaskRequest): Promise<PostTaskResponse> {
    const task = body?.payload?.task

    if (!task) {
      throw new BadRequestException(
        `Expected body.payload.task to be a valid Task, got '${serialize(task)}'`, // prettier-ignore
      )
    }

    const taskId = task?.taskId
    if (!taskId) {
      throw new BadRequestException(
        `Expected body.payload.task.taskId to be a valid ID, got '${serialize(taskId)}'`, // prettier-ignore
      )
    }

    const inputFilenames = await this.fileStoreService.getFilenamesForTask(taskId) // prettier-ignore
    if (!inputFilenames) {
      throw new NotFoundException(`Input files not found for task '${serialize(taskId)}'`)
    }

    await this.taskQueue.emit('tasks', { taskId, inputFilenames }).toPromise()
    return { payload: { taskId } }
  }
}
