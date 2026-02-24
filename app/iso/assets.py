from flask_assets import Bundle

"""
Packages css and js files and applies filters to them like closure compiler
"""

common_css = Bundle(
    'css/vendor/bootstrap.css',
    'css/vendor/bootstrap-grid.css',
    'css/vendor/bootstrap-reboot.css',
    'css/vendor/typeahead.css',
    'css/vendor/helper.css',
    'css/main.css',
    'css/vendor/bootstrap-multiselect.css',
    'css/vendor/prism.css',
    'css/vendor/datatables.css',
    'css/vendor/buttons.dataTables.min.css',
    'css/vendor/font-awesome/css/fontawesome-all.css',
    'css/vendor/jquery.fileupload.css',
    'css/main.css',
    filters='cssmin',
    output='public/css/common.css'
)

common_js = Bundle(
    'js/vendor/jquery.min.js',
    'js/vendor/index.js',
    'js/vendor/popper.min.js',
    'js/vendor/utils.js',
    'js/vendor/tether.js',
    'js/vendor/bootstrap.min.js',
    'js/vendor/util.js',
    'js/vendor/tab.js',
    'js/vendor/button.js',
    'js/vendor/tooltip.js',
    'js/vendor/popover.js',
    'js/vendor/dropdown.js',
    'js/vendor/modal.js',
    'js/vendor/alert.js',
    'js/vendor/bootstrap.bundle.min.js',
    'js/vendor/jquery.validate.min.js',
    'js/vendor/additional-methods.min.js',
    'js/vendor/typeahead.js',
    'js/vendor/bootstrap-multiselect.js',
    'js/vendor/prism.js',
    'js/vendor/datatables.js',
    'js/vendor/highcharts/highcharts.js',
    'js/vendor/highcharts/highcharts-more.js',
    'js/vendor/highcharts/modules/exporting.js',
    'js/vendor/pdfobject.min.js',
    'js/vendor/jvenn.js',
    'js/vendor/canvas2svg.js',
    'js/vendor/dataTables.buttons.min.js',
    'js/vendor/buttons.flash.min.js',
    'js/vendor/jszip.min.js',
    'js/vendor/pdfmake.min.js',
    'js/vendor/vfs_fonts.js',
    'js/vendor/buttons.html5.min.js',
    'js/vendor/buttons.print.min.js',
    'js/vendor/jquery.ui.widget.js',
    'js/vendor/jquery.iframe-transport.js',
    'js/vendor/jquery.fileupload.js',
    Bundle(
        'js/main.js',
        filters='jsmin'
    ),    
    output='public/js/common.js'
)




