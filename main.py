__author__ = 'similarface'
from mthaplogroup_2m_n import main_2m_n
from mthaplogroup_x2 import main_x2
from mthaplogroup_x_point_5 import main_x_point_5
import sys

if __name__ == "__main__":
    filename="/Users/similarface/Documents/FangGE.txt"
    main_2m_n(filename)
    main_x_point_5(filename)
    main_x2(filename)

