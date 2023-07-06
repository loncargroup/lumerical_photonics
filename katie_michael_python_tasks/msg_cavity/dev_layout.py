import numpy as np
import phidl.geometry as pg
from phidl import Device


def make_aperture(x, y, width, height, mask_angle, layer):
    return pg.rectangle((width, height), layer=layer).move(destination=(x - width / 2, y - height / 2)).rotate(
        mask_angle, (x, y))


def make_marker(x0, y0, layer):
    small_sq = 1.75
    large_sq = 1.975

    mark = Device()
    big_rect = pg.rectangle((5, 5), layer=layer).move(destination=(x0 - 2.5, y0 - 2.5))
    mark.add(pg.rectangle((large_sq, large_sq), layer=layer).move(destination=(x0 - 2, y0 - 2)))
    mark.add(pg.rectangle((small_sq, small_sq), layer=layer).move(destination=(x0 - 2, y0 + 2 - small_sq)))
    mark.add(pg.rectangle((large_sq, large_sq), layer=layer).move(destination=(x0 + 2 - large_sq, y0 + 2 - large_sq)))
    mark.add(pg.rectangle((small_sq, small_sq), layer=layer).move(destination=(x0 + 2 - small_sq, y0 - 2)))
    return mark, big_rect


def make_tri(x0, y0, layer):
    tri = pg.taper(20, width1=20, width2=0, layer=layer).rotate(90).move(destination=(x0, y0))
    bigRect = pg.rectangle((25, 25), layer=layer).move(destination=(x0 - 12.5, y0 - 2.5))
    return tri, bigRect


def dev_map(layout):
    layout_index = layout['layout_index']
    aper_index = layout['aper_index']
    layer = layout['layer']
    if layout_index == 0:
        mask_angle = 45
        mask_angle_rad = mask_angle * np.pi / 180
        array_spacer_y = 15
        array_spacer_x = 100
        subarray_spacer_y = 3
        marker_spacer = 420
        subarray_row_num = 5
        mask_size = [.05, .065, .08, .095, .11]

        array_col_num = 5
        array_row_num = 12
        max_row_vals = 36

        marker_x = -130
        marker_y = 10

        pos_list = []
        E = Device()

        for i in range(array_col_num):
            x_pos = i * array_spacer_x
            if mask_angle == 45:  # Set this differently so that we can maximise device packing
                mid_val = (array_col_num - 1) / 2
                rank = abs(i - mid_val)

                array_row_vals = int(-(max_row_vals - array_row_num) / mid_val * rank + max_row_vals)
                y_start = (max_row_vals - array_row_num) / mid_val * (rank - mid_val) * array_spacer_y / 2
            else:
                array_row_vals = array_row_num
                y_start = 0

            for j in range(array_row_vals):  # Make array rows
                y_pos = j * array_spacer_y + y_start
                for k in range(subarray_row_num):  # Make subarray rows
                    # Make the appropriate rotations
                    x = x_pos * np.cos(mask_angle_rad) - (y_pos + k * subarray_spacer_y) * np.sin(mask_angle_rad)
                    y = x_pos * np.sin(mask_angle_rad) + (y_pos + k * subarray_spacer_y) * np.cos(mask_angle_rad)
                    E.add(make_aperture(x, y, mask_size[k], mask_size[k],
                                        mask_angle, layer))  # Actually make the aperture
                    if k == aper_index:
                        pos_list.append((x, y))

        for i in [0, 1]:
            for j in [0, 1]:
                mark, bigRect = make_marker(marker_x + i * marker_spacer, marker_y + j * marker_spacer, layer)
                inverse = pg.boolean(bigRect, mark, 'not',layer = layer)
                E.add(inverse)

        for i in [.35, .65]:
            for j in [0, 1]:
                mark, bigRect = make_marker(marker_x + i * marker_spacer, marker_y + j * (marker_spacer + 55) - 35,
                                            layer)
                inverse = pg.boolean(bigRect, mark, 'not', layer=layer)
                E.add(inverse)

        triangle, bigRect = make_tri(marker_x, marker_y - 45, layer)
        inverse = pg.boolean(bigRect, triangle, 'not', layer=layer)
        E.add(inverse)
        triangle, bigRect = make_tri(marker_x + marker_spacer, marker_y - 45, layer)
        inverse = pg.boolean(bigRect, triangle, 'not', layer=layer)
        E.add(inverse)

        bigAper_x = marker_x + 60
        big_Aper_y = marker_y - 35
        E.add(make_aperture(bigAper_x, big_Aper_y, 10, 10, 0, layer))

        E.add(make_aperture(bigAper_x - 12, big_Aper_y, .06, .06, 0, layer))
        E.add(make_aperture(bigAper_x + 12, big_Aper_y, .06, .06, 0, layer))
        E.add(make_aperture(bigAper_x, big_Aper_y + 12, .06, .06, 0, layer))

        E.add(make_aperture(bigAper_x - 10, big_Aper_y, .1, .1, 0, layer))
        E.add(make_aperture(bigAper_x + 10, big_Aper_y, .1, .1, 0, layer))
        E.add(make_aperture(bigAper_x, big_Aper_y + 10, .1, .1, 0, layer))

        E.add(make_aperture(bigAper_x - 8, big_Aper_y, .2, .2, 0, layer))
        E.add(make_aperture(bigAper_x + 8, big_Aper_y, .2, .2, 0, layer))
        E.add(make_aperture(bigAper_x, big_Aper_y + 8, .2, .2, 0, layer))

        x_shift = -marker_x - 210
        y_shift = -marker_y - 200
        E.move(origin=(x_shift, y_shift))
        pos_list = [(i + x_shift, j + y_shift) for (i, j) in pos_list]

    return E, pos_list
