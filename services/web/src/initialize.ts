import type { Router } from 'next/router'

import { configureStore } from 'src/state/store'
import { getTaskIdTrigger } from 'src/state/task/task.actions'

export interface InitializeParams {
  router: Router
}

export async function initialize({ router }: InitializeParams) {
  void router.prefetch('/') // eslint-disable-line no-void
  void router.prefetch('/results') // eslint-disable-line no-void

  const { persistor, store } = await configureStore({ router })

  store.dispatch(getTaskIdTrigger())

  return { persistor, store }
}
