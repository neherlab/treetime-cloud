import { Injectable } from '@nestjs/common'

import uuidv4 from 'uuid/v4'

@Injectable()
export class TaskService {
  public generateTaskId() {
    return uuidv4()
  }
}
