import actionCreatorFactory from 'typescript-fsa'

const action = actionCreatorFactory('WS')

/* Action creator: websocket saga failed with an error */
export const wsFailed = action<{ error: string }>('FAILED')

/* Action creator: websocket saga was cancelled */
export const wsCancelled = action<void>('CANCELLED')
