import React from 'react'

import FileUploadZone from '../components/FileUpload/FileUploadZone'

import './Home.scss'

function Home() {
  return (
    <>
      <h1 className="h1-home">{'Upload files'}</h1>

      {/* prettier-ignore */}
      <div>
        <FileUploadZone/>
      </div>
    </>
  )
}

export default Home
