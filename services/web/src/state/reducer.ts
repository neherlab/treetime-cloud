import { combineReducers } from 'redux'
import { persistReducer } from 'redux-persist'
import storage from 'redux-persist/lib/storage'
import { routerReducer } from 'connected-next-router'
import { RouterState } from 'connected-next-router/types'

import { TaskState } from './task/task.state'
import { taskReducer } from './task/task.reducer'

import { SettingsState } from './settings/settings.state'
import { settingsReducer } from './settings/settings.reducer'

import { UploadState } from './upload/upload.state'
import { uploadReducer } from './upload/upload.reducer'

export interface State {
  settings: SettingsState
  task: TaskState
  router: RouterState
  upload: UploadState
}

const SETTINGS_VERSION = 1
const settingsReducerPersisted = persistReducer(
  { key: 'settings', version: SETTINGS_VERSION, storage, timeout: 3000 },
  settingsReducer,
)

const rootReducer = () =>
  combineReducers({
    router: routerReducer,
    settings: settingsReducerPersisted,
    task: taskReducer,
    upload: uploadReducer,
  })

export default rootReducer
