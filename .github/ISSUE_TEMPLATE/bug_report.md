---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: ''
assignees: ''

---

**The code version and commit ID that you are using**
- (How to find? In your local code, run the following command, and copy the output here)
- `echo $(git rev-parse --abbrev-ref HEAD) $(git rev-parse --short HEAD)`

**Describe the bug**
- What is the bug?
- What you expected to happen without the bug?

**Your `config.log` file**
- The `config.log` file in your code root directory can help developers solving most simple problems. Unless you are sure the compilation environment will not help developers, you are encouraged to attach it.
- **Otherwise**, please mentioned your platform and environment: 
  - OS and version: [e.g. Ubuntu 24]
  - Compiler and version of Fortran, C++, MPI and CUDA (if applicable)
