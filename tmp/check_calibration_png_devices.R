source('tmp/preview_current_calibration.R')
for (device_type in c('cairo-png', 'windows')) {
  png(file.path(e$output_dir, paste0('device_', device_type, '.png')),
      width = 16, height = 11.4, units = 'in', res = 180, type = device_type)
  e$draw_figure('pip_calibration', e$scenario_groups$pure, methods = 'SuSiE-slide')
  dev.off()
}
