This is a repository of the `prog/naoUser` folder in the lab files initially given to us. The important file is `src/balance_task.cpp`. You can either remove your existing naoUser folder and replace it with this one by doing this:
```bash
cd prog/
rm -rf naoUser
git clone https://github.com/alexanderkoumis/545-Project.git naoUser
```

Or, change the remote of your existing naoUser folder by doing this:
```bash
cd prog/naoUser
git stash
git remote set-url origin https://github.com/alexanderkoumis/545-Project.git
git pull
```

Either way, you will need to make the build folder again and run cmake like you did when setting up the simulator:
```
BUILD_DIR=naoUser/x86_64 # or naoUser/x86_64mac if you're on mac
cd $BUILD_DIR
cmake ../src
cd ..
make install
```

Another option is to just copy `src/balance_task.cpp` to your folder without cloning the project. You would then have to modify `src/initUserTasks.c` and `src/CMakeLists.txt` accordingly. Then, when you want to commit your changes to the repo, you would clone this repo and copy over your balance_task.cpp, commit it, and push.