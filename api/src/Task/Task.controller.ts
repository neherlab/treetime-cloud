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

import { UploadUnique } from '../common/UploadUnique'

import { FileStoreService, UploadedFileData } from './FileStore.service'
import { TaskIdService } from './TaskId.service'

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
  public async postTask(
    @Body() { payload: { task } }: PostTaskRequest,
  ): Promise<PostTaskResponse> {
    const taskId = task.id
    const inputFilepaths = await this.fileStoreService.getFilepathsForTask(taskId) // prettier-ignore
    await this.taskQueue.emit('tasks', { taskId, inputFilepaths }).toPromise()
    return { payload: { taskId } }
  }
}
