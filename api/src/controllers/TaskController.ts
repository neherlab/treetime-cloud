import { Body, JsonController, Post, UploadedFile } from 'routing-controllers'

import { File } from 'multer'

export interface Task {
  id: string
}

export interface UploadRequestBody {
  taskId: string
}

export interface PostTaskRequest {
  payload: {
    task: Task
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
}
