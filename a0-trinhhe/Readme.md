# Assignment 0 - Settings things up

This assignment has two purposes:

1. match your github username to your student ID
2. make sure you can build and run the starter code

## Part 1: match github to student ID

1. Insert your information (full name, student id, nethz username) into `my-info.json`.

2. Commit and push the changes to `my-info.json`. For example:

   ```
   git add my-info.json
   git commit -m "added my information"
   git push
   ```

3. Verify that your changes were uploaded to github.com by checking github.com/cmm-21/a0-YOUR_GITHUB_NAME

## Part 2: build and run starter code

This assignement is to make sure your computer and you are both set up for future assignements. For future assignements, you will have to commit and push your solutions to github. Commits (or rather pushes) after the deadline will be disregarded. So make sure you are ready and test out things with this assignment.

### Prerequisites
Make sure you install the following:

- Git (https://git-scm.com)
    + Windows: download installer from website
    + Linux: e.g. `sudo apt install git`
    + MacOS: e.g. `brew install git`
- CMake (https://cmake.org/) <= 3.19
    + Windows: download installer from website
    + Linux: e.g. `sudo apt install cmake`
    + MacOS: e.g. `brew install cmake`
    + **Note** the latest CMake (3.20) is not compatible with our skeleton code: see issue [i4](https://github.com/cmm-21/a0/issues/4) and [i5](https://github.com/cmm-21/a0/issues/5).

### Building the code

On Windows, you can use Git Bash to perform the steps mentioned below.

1. Git clone your project, e.g.:

   ```
   git clone git@github.com:cmm-21/a0-XXX.git
   ```

2. Run CMake to create project files
    ```
    cd a0-XXX
    mkdir build && cd build
    cmake ..
    ```
3. Compile code and run executable
    - for MacOS and Linux:
        ```
        make
        ./src/app/app
        ```
    - for Windows: 
        * open `build/a0.sln` in Visual Studio
        * in the project explorer, right-click target "app" and set as startup app.
        * Hit F5!

3. Implement the function `math::add` in the file `src/math/add.cpp`.

4. Check if test passed: `make && ctest -V`

5. Commit and push changes, and check on github if test also passed.

### Some comments

- If you are new to git, there are many tutorials online, e.g. http://rogerdudler.github.io/git-guide/.

- You will not have to edit any CMake file, so no real understanding of CMake is required. You might however want to generate build files for your favorite editor/IDE: https://cmake.org/cmake/help/latest/manual/cmake-generators.7.html

- If you want to use [VS Code](https://code.visualstudio.com/), follow the instructions below:
    - Install the following VS Code extensions.
        - [C/C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools)
        - [CMake Tools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools)
    - Select a kit/variant and configure the project. See [this](https://vector-of-bool.github.io/docs/vscode-cmake-tools/getting_started.html#configuring-your-project) for more details (this corresponds to ```cmake ..```)
    - Click "Build" button in status bar (bar at the botton) to compile and run the executable by "Launch" button.
    - If you want to debug your code, set cmake mode into "CMake: [Debug]" and click "Launch the debugger" button. 
        - In "CMake: [Debug]" mode, the executable tends to run very slowly. Select "CMake: [Release]" to run the code fast (but in Release mode you cannot use breakpoints). 

![status bar](vscode.png)