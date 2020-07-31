import actionCreatorFactory from 'typescript-fsa'

const action = actionCreatorFactory('Error')

export interface GenericErrorParams {
  error: Error
}

export const errorAdd = action<GenericErrorParams>('errorAdd')

export const errorDismiss = action<unknown>('errorDismiss')
