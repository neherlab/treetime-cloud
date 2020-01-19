import MiniCssExtractPlugin from 'mini-css-extract-plugin'
import sass from 'sass'

import { findModuleRoot } from '../../../lib/findModuleRoot'

const { moduleRoot } = findModuleRoot()

export interface WebpackLoadStylesParams {
  isDev: boolean
  isServer: boolean
  sourceMap: boolean
  exclude?: string[] | RegExp | RegExp[]
  include?: string[] | RegExp | RegExp[]
  modules?: boolean
}

export default function webpackLoadStyles({
  isDev,
  isServer,
  sourceMap,
  exclude,
  include,
  modules,
}: WebpackLoadStylesParams) {
  const isDevWeb = isDev && !isServer

  return [
    {
      test: /\.(css|sass|scss)$/,
      exclude,
      include,
      sideEffects: true,
      use: [
        isDevWeb && { loader: 'css-hot-loader' },
        { loader: MiniCssExtractPlugin.loader },
        { loader: 'css-loader', options: { modules, sourceMap } },
        {
          loader: 'postcss-loader',
          options: { sourceMap, config: { path: moduleRoot } },
        },
        {
          loader: 'sass-loader',
          options: {
            implementation: sass,
            sourceMap,
            sassOptions: {
              sourceMap,
            },
          },
        },
      ].filter(Boolean),
    },
  ]
}
