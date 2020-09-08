import { EventChannel, eventChannel } from 'redux-saga'
import { call, cancelled, put, take } from 'redux-saga/effects'
import { ActionCreator } from 'typescript-fsa'

import { wsCancelled, wsFailed } from '../state/websockets/websockets.actions'

type WsMessageType = string

/**
 * Creates a redux-saga channel that is subscribed to given message types on
 * a given sockets.io socket.
 *
 * Incoming messages result in calling channel's `emit()`
 * function . Saga can then retrieve the message from the
 * channel by calling `take()` on it.
 *
 */
function createWsChannel(socket: SocketIOClient.Socket, messageTypes: WsMessageType[]) {
  return eventChannel((emit) => {
    // For every message type, subscribe channel's emit to the socket.
    messageTypes.forEach((messageType) => {
      socket.on(messageType, ({ payload }: { payload: unknown }) => {
        emit({ messageType, payload })
      })
    })

    // Cleanup handler: for every message type unsubscribe from the socket.
    return () => {
      messageTypes.forEach((messageType) => {
        socket.off(messageType)
      })
    }
  })
}

/**
 * Creates a saga that, given a sockets.io connection and an array of
 * FSA-compliant redux action creators, listens for incoming
 * socket.io messages, converts them into appropriate redux actions and
 * dispatches these actions to redux store.
 */
export default function createWsReceiverSaga(
  socket: SocketIOClient.Socket,
  actionCreators: ActionCreator<any>[], // eslint-disable-line @typescript-eslint/no-explicit-any
) {
  const actionTypes = actionCreators.map((ac) => ac.type)

  return function* wsReceiverSaga() {
    const channel = (yield call(createWsChannel, socket, actionTypes)) as EventChannel<unknown>

    while (true) {
      try {
        const { messageType, payload } = (yield take(channel)) as { messageType: WsMessageType; payload: unknown }
        const actionCreator = actionCreators.find((ac) => ac.type === messageType)
        if (actionCreator) {
          yield put(actionCreator(payload))
        }
      } catch (error) {
        yield put(wsFailed({ error }))
        throw error
      } finally {
        if (yield cancelled()) {
          channel.close()
          yield put(wsCancelled())
        }
      }
    }
  }
}
