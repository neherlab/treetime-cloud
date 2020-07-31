import React, { PropsWithChildren, HTMLProps } from 'react'

import styled from 'styled-components'

import { NavigationBar } from './NavigationBar'

export const Container = styled.div`
  max-width: 1700px;
  margin: 0 auto;
`

const Header = styled.header``

const MainContent = styled.main`
  margin: 0 1rem;
`

export type LayoutMainProps = PropsWithChildren<HTMLProps<HTMLDivElement>>

export function Layout({ children }: LayoutMainProps) {
  return (
    <Container>
      <Header>
        <NavigationBar />
      </Header>

      <MainContent>{children}</MainContent>
    </Container>
  )
}
