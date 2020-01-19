const path = require('path')

module.exports = {
  root: true,
  parser: '@typescript-eslint/parser',
  parserOptions: {
    sourceType: 'module',
    ecmaFeatures: {
      jsx: false,
      globalReturn: false,
    },
    project: './tsconfig.json',
    extraFileExtensions: ['.json'],
    warnOnUnsupportedTypeScriptVersion: true,
  },
  globals: {},
  extends: [
    'eslint:recommended',
    'airbnb-base',

    'plugin:@typescript-eslint/eslint-recommended',
    'plugin:@typescript-eslint/recommended',
    'plugin:@typescript-eslint/recommended-requiring-type-checking',
    'plugin:array-func/all',
    'plugin:import/errors',
    'plugin:import/typescript',
    'plugin:import/warnings',
    'plugin:jest/recommended',
    'plugin:jest/style',
    'plugin:lodash/recommended',
    'plugin:node/recommended',
    'plugin:promise/recommended',
    'plugin:security/recommended',
    'plugin:sonarjs/recommended',
    'plugin:unicorn/recommended',

    // prettier should go last
    'plugin:prettier/recommended',
    'prettier/@typescript-eslint',
  ],
  plugins: [
    'array-func',
    'cflint',
    'import',
    'jest',
    'json',
    'lodash',
    'no-loops',
    'no-secrets',
    'node',
    'only-ascii',
    'promise',
    'security',
    'sonarjs',
    'unicorn',

    'only-warn',

    '@typescript-eslint',
    '@typescript-eslint/tslint',

    // prettier should go last
    'prettier',
  ],
  rules: {
    '@typescript-eslint/array-type': 'off',
    '@typescript-eslint/explicit-function-return-type': 'off',
    'array-func/prefer-array-from': 'off',
    'cflint/no-substr': 'warn',
    'cflint/no-this-assignment': 'warn',
    'import/extensions': [
      'warn',
      'ignorePackages',
      {
        js: 'never',
        jsx: 'never',
        mjs: 'never',
        ts: 'never',
        tsx: 'never',
      },
    ],
    'import/no-extraneous-dependencies': ['warn', { devDependencies: true }],
    'import/order': 'warn',
    'import/prefer-default-export': 'off',
    'jest/consistent-test-it': 'warn',
    'jest/expect-expect': 'warn',
    'jest/no-test-callback': 'warn',
    'lodash/chaining': 'off',
    'lodash/import-scope': 'off',
    'lodash/prefer-constant': 'off',
    'lodash/prefer-lodash-chain': 'off',
    'lodash/prefer-lodash-method': 'off',
    'no-console': ['warn', { allow: ['info', 'warn', 'error'] }],
    'no-loops/no-loops': 'warn',
    'no-secrets/no-secrets': ['warn', { tolerance: 5 }],
    'node/no-missing-import': 'off',
    'node/no-unpublished-import': 'off',
    'node/no-unpublished-require': 'off',
    'node/no-unsupported-features/es-builtins': 'off',
    'node/no-unsupported-features/es-syntax': 'off',
    'only-ascii/only-ascii': 'warn',
    'prettier/prettier': 'warn',
    'security/detect-non-literal-fs-filename': 'off',
    'security/detect-object-injection': 'off',
    'unicorn/filename-case': 'off',
    'unicorn/new-for-builtins': 'off',
    'unicorn/no-abusive-eslint-disable': 'warn',
    'unicorn/prefer-query-selector': 'off',
    'unicorn/prevent-abbreviations': 'off',

    'lines-between-class-members': [
      'warn',
      'always',
      { exceptAfterSingleLine: true },
    ],

    '@typescript-eslint/tslint/config': [
      'warn',
      { lintFile: path.join(__dirname, 'tslint.json') },
    ],

    'require-await': 'off',
    '@typescript-eslint/require-await': 'off',

    'no-unused-expressions': 'off',
    '@typescript-eslint/no-unused-expressions': 'warn',
  },
  env: {
    browser: false,
    es6: true,
    jest: true,
    node: true,
  },
  overrides: [
    {
      files: ['*.json'],
      rules: {
        '@typescript-eslint/no-useless-files': 'off',
      },
    },
    {
      files: ['*.d.ts'],
      rules: {
        '@typescript-eslint/no-unused-vars': 'off',
        'no-useless-constructor': 'off',
      },
    },
    {
      files: [
        '.eslintrc.js',
        'babel.config.js',
        'config/**/*.js',
        'config/**/*.ts',
      ],
      rules: {
        '@typescript-eslint/no-var-requires': 'off',
        'sonarjs/cognitive-complexity': ['warn', 50],
      },
    },
    {
      files: ['**/src/**/__tests__/integration/*.test.*'],
      rules: {
        'jest/expect-expect': 'off',
      },
    },
    {
      files: ['**/src/controllers/**'],
      rules: {
        'class-methods-use-this': 'off',
      },
    },
  ],
}
