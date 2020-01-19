import { JsonController, Post } from 'routing-controllers'

@JsonController()
export default class TaskController {
  @Post('/api/v1/task')
  public async postTask() {
    return { payload: {} }
  }
}
