import {
  Body,
  Get,
  JsonController,
  Post,
  Req,
  UseBefore,
} from 'routing-controllers'

import uuidv4 from 'uuid/v4'

import multer from 'multer'

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

export interface RequestWithFiles {
  files: {
    DATES: File[]
    FASTA: File[]
    NWK: File[]
  }
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

  @Post('/api/v1/upload')
  @UseBefore(
    multer().fields([
      { name: 'DATES', maxCount: 1 },
      { name: 'FASTA', maxCount: 1 },
      { name: 'NWK', maxCount: 1 },
    ]),
  )
  public async uploadFiles(
    @Body({ required: true }) { taskId }: UploadRequestBody,
    @Req() req: RequestWithFiles,
  ): Promise<PostTaskResponse> {
    // TODO: prevent crash when a file of unknown type is uploaded

    const dates = req?.files?.DATES
    const fasta = req?.files?.FASTA
    const nwk = req?.files?.NWK

    if (dates) {
      this.files.set('DATES', dates[0])
    }

    if (fasta) {
      this.files.set('FASTA', fasta[0])
    }

    if (nwk) {
      this.files.set('NWK', nwk[0])
    }

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
