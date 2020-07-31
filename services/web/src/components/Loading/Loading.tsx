import React from 'react'

import { useTranslation } from 'react-i18next'

import { ReactComponent as Logo } from 'src/assets/img/logo.svg'
import styled from 'styled-components'

const Container = styled.div`
  display: flex;
  width: 100%;
  height: 100%;
  overflow: hidden;
`

const SpinningLogo = styled(Logo)`
  margin: auto;
  width: 80px;
  height: 80px;
  animation: spin 0.5s linear infinite;
  @keyframes spin {
    100% {
      transform: rotate(360deg);
    }
  }
`

function Loading() {
  const { t } = useTranslation()
  return (
    <Container title={t('Loading...')}>
      <SpinningLogo />
    </Container>
  )
}

export default Loading
