import { Injectable } from '@nestjs/common'

import uuidv4 from 'uuid/v4'

@Injectable()
export class TaskIdService {
  public generateTaskId() {
    return uuidv4()
  }
}
