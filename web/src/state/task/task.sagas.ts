import { takeLatest } from 'redux-saga/effects'

import axios from 'axios'

import { getTaskIdAsync, triggerGetTaskId } from './task.actions'

import fsaSaga from '../util/fsaSaga'

function getTaskApi() {
  return axios.get(`/api/v1/taskId`, {
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'multipart/form-data',
    },
  })
}

export const taskSaga = takeLatest(
  triggerGetTaskId,
  fsaSaga(getTaskIdAsync, getTaskApi),
)

export default [taskSaga]
