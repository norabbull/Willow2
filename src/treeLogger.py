# # -*- coding: utf-8 -*-
# """
# Created on Mon Jun 28 12:26:02 2021

# @author: norab
# """

# import logging
# # from datetime import datetime
# # import os

# FORMATTER = logging.Formatter("%(asctime)s — %(name)s — %(levelname)s — %(message)s")
# LOGFILE = 'Willow.log'
# # log_dir = os.path.join(os.path.normpath(os.getcwd() + os.sep + os.pardir), 'logs')
# # LOGFILE = os.path.join(log_dir, 'Willow.log')


# def get_console_handler():
#     console_handler = logging.StreamHandler()
#     console_handler.setFormatter(FORMATTER)
#     return console_handler


# def get_file_handler():
#     file_handler = logging.FileHandler(filename=LOGFILE)
#     file_handler.setFormatter(FORMATTER)
#     return file_handler


# def get_logger(logger_name):
    
#     # if 'datetime' in str(path): 
#     #     path = path.replace('datetime', datetime.now().strftime("%d.%m.%Y_%H.%M"))

#     logger = logging.getLogger(logger_name)
#     logger.setLevel(level=logging.DEBUG) # better to have too much log than not enough
#     logger.addHandler(get_file_handler())
#     logger.addHandler(get_console_handler())
    
#     # with this pattern, it's rarely necessary to propagate the error up to parent
#     logger.propagate = False
    
#     return logger











#%%

logging.basicConfig(level=logging.DEBUG, 
                     format='%(asctime)s - %(name)s - %(levelname)s - %message)s',
                     datefmt='%m/%d/%Y %H:%M:&S')

logger = logging.getLogger(__name__)

