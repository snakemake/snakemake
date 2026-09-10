#!/usr/bin/env bash
#
set -exu -o pipefail

echo -e "#!/usr/bin/env bash\necho fake script installed!" >fake
chmod +x fake
mv fake "${CONDA_PREFIX}"/bin/.
