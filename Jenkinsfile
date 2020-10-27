//----------------------------------------------------------------------------------------------------------------------
// This declarative Jenkins pipeline encodes all the steps required for the nightly/continuous of a single platform.
// Other jobs may call this pipeline to execute the build, test and installation of a set platforms.
//
// Author: Pere Mato
//----------------------------------------------------------------------------------------------------------------------

pipeline {
  parameters {
    string(name: 'EXTERNALS', defaultValue: 'devgeantv/latest', description: 'LCG software stack in CVMFS')
    choice(name: 'MODE', choices: ['experimental', 'nightly', 'continuous'], description: 'CDash mode')
    string(name: 'ExtraCMakeOptions', defaultValue: '', description: 'CMake extra configuration options')
    string(name: 'LABEL', defaultValue: 'centos7', description: 'Jenkins label for physical nodes or container image for docker')
    choice(name: 'COMPILER', choices: ['gcc7', 'gcc8', 'gcc9', 'gcc10', 'clang8', 'clang10', 'native'])
    choice(name: 'BUILDTYPE', choices: ['Release', 'Debug'])
    choice(name: 'OPTION', choices: ['default', 'SPEC', 'AVX', 'GDML'])
    choice(name: 'BACKEND', choices: ['scalar', 'vc'])
    string(name: 'DOCKER_LABEL', defaultValue: 'docker-host-noafs', description: 'Label for the the nodes able to launch docker images')
  }

  environment {
    CMAKE_INSTALL_PREFIX = 'install'
    CMAKE_SOURCE_DIR     = 'vecgeom'
    CMAKE_BINARY_DIR     = 'build'
  }

  agent none

  stages {
    //------------------------------------------------------------------------------------------------------------------
    //---Build & Test stages--------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------------------
    stage('Prepa'){
      steps {
        init()
      }
    }
    stage('InDocker') {
      when {
        beforeAgent true
        expression { params.LABEL =~ 'centos|ubuntu' && !(params.LABEL =~ 'physical')}
      }
      agent {
        docker {
          image "gitlab-registry.cern.ch/sft/docker/$LABEL"
          label "$DOCKER_LABEL"
          args  """-v /cvmfs:/cvmfs 
                   -v /ccache:/ccache 
                   -v /ec:/ec
                   -e SHELL -e gitlabMergedByUser -e gitlabMergeRequestIid
                   --net=host
                   --hostname ${LABEL}-docker
                """
        }
      }
      stages {
        stage('Build&Test') {
          steps {
            buildAndTest()
          }
          post {
            success {
              deleteDir()
            }
          }
        }
      }
    }
    stage('InGPU') {
      when {
        beforeAgent true
        expression { params.LABEL =~ 'cuda|physical' }
      }
      agent {
        label 'cuda10'
      }
      stages {
        stage('Build&Test') {
          steps {
            buildAndTest()
          }
          post {
            success {
              deleteDir()
            }
          }
        }
      }
    }
  }
}

def init() {
  currentBuild.displayName = "#${BUILD_NUMBER}" + ' ' + params.OPTION + '-' + params.BACKEND + '-' + params.LABEL + '-' + params.COMPILER + '-' + params.BUILDTYPE
}

def buildAndTest() {
  sh label: 'build_and_test', script: """
    source /cvmfs/sft.cern.ch/lcg/views/${EXTERNALS}/x86_64-centos7-${COMPILER}-opt/setup.sh
    env | sort | sed 's/:/:?     /g' | tr '?' '\n'
    ctest -VV -S vecgeom/jenkins/vecgeom-ctest.cmake,$MODE
  """
}