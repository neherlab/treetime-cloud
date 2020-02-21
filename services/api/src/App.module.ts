import { Module } from '@nestjs/common'

import { TaskModule } from './Task/Task.module'

@Module({
  imports: [TaskModule],
})
export class AppModule {}
