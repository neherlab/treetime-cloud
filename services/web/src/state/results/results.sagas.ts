import _ from 'lodash'

import { fork } from 'redux-saga/effects'
import io from 'socket.io-client'

import createWsReceiverSaga from 'src/state/util/websocketReceiverSaga'

const host = process.env.BACKEND_HOST_EXTERNAL
const port = process.env.BACKEND_PORT_EXTERNAL
if (!host || !port) {
  throw new Error(`Agent host/port invalid: host: '${host}', port: ${port}`)
}

const url = `${host}:${port}`
const socket = io(url, { path: '/websockets' })
const actions = _.uniq([getFindingsForSeriesSuccessful, getFindingsForSeriesFailed])

export default [fork(createWsReceiverSaga(socket, actions))]
