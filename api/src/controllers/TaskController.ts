import {
  Body,
  Get,
  JsonController,
  Post,
  UploadedFile,
} from 'routing-controllers'

import uuidv4 from 'uuid/v4'

import { File } from 'multer'

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

@JsonController()
export default class TaskController {
  // HACK: should become a service, with client isolation
  private files: Map<string, File> = new Map<string, File>()

  @Get('/api/v1/taskId')
  public async getTaskId(): Promise<GetTaskIdResponse> {
    const taskId = this.generateTaskId()
    return { payload: { taskId } }
  }

  @Post('/api/v1/upload/dates')
  public async uploadDates(
    @Body({ required: true }) { taskId }: UploadRequestBody,
    @UploadedFile('file', { required: true }) dates: File,
  ): Promise<PostTaskResponse> {
    this.files.set('DATES', dates)
    return { payload: { taskId } }
  }

  @Post('/api/v1/upload/fasta')
  public async uploadFasta(
    @Body({ required: true }) { taskId }: UploadRequestBody,
    @UploadedFile('file', { required: true }) fasta: File,
  ): Promise<PostTaskResponse> {
    this.files.set('FASTA', fasta)
    return { payload: { taskId } }
  }

  @Post('/api/v1/upload/nwk')
  public async uploadNwk(
    @Body({ required: true }) { taskId }: UploadRequestBody,
    @UploadedFile('file', { required: true }) nwk: File,
  ): Promise<PostTaskResponse> {
    this.files.set('NWK', nwk)
    return { payload: { taskId } }
  }

  @Post('/api/v1/task')
  public async postTask(
    @Body({ required: true })
    { payload: { task } }: PostTaskRequest,
  ): Promise<PostTaskResponse> {
    return { payload: { taskId: task.id } }
  }

  // TODO: should become a part of a new service
  private generateTaskId() {
    return uuidv4()
  }
}
