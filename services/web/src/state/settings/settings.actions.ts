import actionCreatorFactory from 'typescript-fsa'

import type { LocaleKey } from 'src/i18n/i18n'

const action = actionCreatorFactory('Settings')

export const setLocale = action<LocaleKey>('setLocale')
