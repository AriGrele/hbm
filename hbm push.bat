@echo off
set /P mes=Message? 
git add .
git commit -m %mes%
git push origin HEAD:master