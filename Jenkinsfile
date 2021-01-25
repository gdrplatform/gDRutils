DockerfilePipeline(
  preBuild: { tags,labels, rpServiceUrl, rpCookie ->
    def rpsession_cookie = "${rpCookie}"
    def rp_service_url = "${rpServiceUrl}"
    def rp_project_id = "604"
    def GIT_URL = env.GIT_URL ? env.GIT_URL : env.GIT_URL_1

    // create new tag from branch's name (and commit for devel and master)
    def tag = null
    def commit_11 = sh returnStdout: true, script: "echo ${GIT_COMMIT} | head -c 11"
    if ("${GIT_BRANCH}" == "master") {
        tag = "${GIT_BRANCH}-${commit_11}"
        tag = tag.trim()
    } else {
        if ("${GIT_BRANCH}" == "devel") {
            tag = "${GIT_BRANCH}-${commit_11}"
        } else {
            tag = sh returnStdout: true, script: "sed 's=[^[:alnum:]-]=-=g' <<< ${GIT_BRANCH}"
        }
        tag = tag.trim()
        labels.registry_credentials = "rplatform_prodreg-credentials"
        labels.temp_registry = "registry.rplatform.org:5000"
    }
    tags.add(tag)
    echo "Tag docker image with ${tag}"
    response = sh(script: """
        curl --retry 5 -w '%{http_code}' --retry-delay 0 -H "Cookie: rpsession=${rpsession_cookie}" -H "Content-Type: application/json" --data '{"commit":"${GIT_COMMIT}", "gitRepositoryUrl":"${GIT_URL}", "tag":"${tag}", "branch":"${GIT_BRANCH}", "useCache":"false", "useDepth":"0", "cloneSubmodules":"false"}' -X POST "${rp_service_url}/v1/build"
    """, returnStdout: true)
    if ( "${response}".contains("200") ) {
        echo "RP service was notified that the build is started."
    } else {
        echo "Failing build. Cause: Something went wrong with API request (timeout, bad response). Error code: ${response}"
        error("Failing build. Cause: Something went wrong with API request (timeout, bad response)")
    }

    if (!fileExists("Dockerfile")) {
        echo "Moving ./rplatform/Dockerfile to ./Dockerfile"
        sh 'cp ./rplatform/Dockerfile ./Dockerfile'
        sh 'ls -la'
    }

  },
  onFailure: { tags,labels, rpServiceUrl, rpCookie ->
    def rpsession_cookie = "${rpCookie}"
    def rp_service_url = "${rpServiceUrl}"
    def rp_project_id = "604"

    sh 'export PYTHONWARNINGS="ignore:Unverified HTTPS request"'
    sh 'printenv | sort'
    response = sh(script: """
        curl --retry 5 -w '%{http_code}' --retry-delay 0 -X PUT -H "Cookie: rpsession=${rpsession_cookie}" "${rp_service_url}/v1/statuses?status=failed&commit=${GIT_COMMIT}"
    """, returnStdout: true)
    if ( "${response}".contains("200") ) {
        echo "Image status updated."
    } else {
        echo "Failing build. Cause: Something went wrong with API request (timeout, bad response). Error code: ${response}"
        error("Failing build. Cause: Something went wrong with API request (timeout, bad response)")
    }
  },
  onSuccess: { tags,labels,rpServiceUrl,rpCookie ->
    def rpsession_cookie = "${rpCookie}"
    def rp_service_url = "${rpServiceUrl}"
    def rp_project_id = "604"
    sh 'printenv | sort'

    // overwrite docker credentials and registry address
    if ("${GIT_BRANCH}" != "master") {
        labels.registry_credentials = "rplatform_prodreg-credentials"
        labels.temp_registry = "registry.rplatform.org:5000"
        registry = "registry.rplatform.org:5000"
    }
    echo "Registry: ${registry}"

    // create new tag from branch's name (and commit for devel and master)
    def tag = null
    def image = null
    def cbs_docker_image = "${registry}/githubroche/gdrplatform/gdrutils"
    def commit_11 = sh returnStdout: true, script: "echo ${GIT_COMMIT} | head -c 11"
    if ("${GIT_BRANCH}" == "master" || "${GIT_BRANCH}" == "devel") {
        // master goes to CBS registry
        tag = "${GIT_BRANCH}-${commit_11}"
    } else {
        // devel and other branches go to RP temp registry
        if ("${GIT_BRANCH}" == "devel") {
            tag = "${GIT_BRANCH}-${commit_11}"
        } else {
            tag = sh returnStdout: true, script: "sed 's=[^[:alnum:]-]=-=g' <<< ${GIT_BRANCH}"
        }
        tag = tag.trim()
        echo "Tag docker image with ${tag}"
        docker.withRegistry("https://${registry}", labels.registry_credentials) {
            sh "docker tag ${cbs_docker_image}:${env.BUILD_ID} ${cbs_docker_image}:${tag}"
            sh "docker tag ${cbs_docker_image}:${env.BUILD_ID} ${cbs_docker_image}:latest"
            echo "Pushing ${tag}"
            image = docker.image("${cbs_docker_image}:${tag}")
            image.push()
            echo "Pushing latest"
            image = docker.image("${cbs_docker_image}:latest")
            image.push()
        }
    }

    sh 'docker images'
    sh "ls -la ${WORKSPACE}"

    // notify service about success build
    response = sh(script: """
        curl --retry 5 -w '%{http_code}' --retry-delay 0 -X PUT -H "Cookie: rpsession=${rpsession_cookie}" "${rp_service_url}/v1/statuses?status=success&commit=${GIT_COMMIT}"
    """, returnStdout: true)
    if ( "${response}".contains("200") ) {
        echo "Image status updated."
    } else {
        echo "Failing build. Cause: Something went wrong with API request (timeout, bad response). Error code: ${response}"
        error("Failing build. Cause: Something went wrong with API request (timeout, bad response)")
    }

    // run tests and notify rp service
    try {
        // tests IN_PROGRESS
        response = sh(script: """
          curl --retry 5 -w '%{http_code}' --retry-delay 0 -H 'Cookie: rpsession=${rpsession_cookie}' \
          -X PUT "${rp_service_url}/v1/statuses?status=in_progress&branch_name=${GIT_BRANCH}&project_id=${rp_project_id}&build_log=${BUILD_URL}"
        """, returnStdout: true)
        if ( "${response}".contains("200") ) {
          echo "Image status updated."
        } else {
          echo "Failing build. Cause: Something went wrong with API request (timeout, bad response). Error code: ${response}"
          error("Failing build. Cause: Something went wrong with API request (timeout, bad response)")
        }
        // docker RUN TESTS
        sh "docker run --rm -v ${WORKSPACE}:/mnt/vol -i -e PASSWORD=run_tests -w='/mnt/vol' ${cbs_docker_image}:${env.BUILD_ID} /bin/sh -c 'rp_env=DOCKER /bin/sh rplatform/run_tests.sh'"
        // tests SUCCESS
        response = sh(script: """
          curl --retry 5 -w '%{http_code}' --retry-delay 0 -H "Cookie: rpsession=${rpsession_cookie}" \
            -X PUT "${rp_service_url}/v1/statuses?status=success&branch_name=${GIT_BRANCH}&project_id=${rp_project_id}&build_log=${BUILD_URL}"
        """, returnStdout: true)
        if ( "${response}".contains("200") ) {
            echo "Image status updated."
        } else {
            echo "Failing build. Cause: Something went wrong with API request (timeout, bad response). Error code: ${response}"
            error("Failing build. Cause: Something went wrong with API request (timeout, bad response)")
        }
    } catch (err) {
        echo "Notifying repo with FAIL"
        response = sh(script: """
          curl --retry 5 -w '%{http_code}' --retry-delay 0 -H "Cookie: rpsession=${rpsession_cookie}" \
              -X PUT "${rp_service_url}/v1/statuses?status=failed&branch_name=${GIT_BRANCH}&project_id=${rp_project_id}&build_log=${BUILD_URL}"
        """, returnStdout: true)
        if ( "${response}".contains("200") ) {
            echo "Image status updated."
        } else {
            echo "Failing build. Cause: Something went wrong with API request (timeout, bad response). Error code: ${response}"
            error("Failing build. Cause: Something went wrong with API request (timeout, bad response)")
        }
        throw new Exception("Tests failed: ${err}.")
    }
  }
)