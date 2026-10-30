@echo off

SET BASE_PATH=%~dp0

CALL julia +1.11.9 --project=%BASE_PATH% %BASE_PATH%\format.jl