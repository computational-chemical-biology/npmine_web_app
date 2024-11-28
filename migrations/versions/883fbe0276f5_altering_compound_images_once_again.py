"""altering compound_images once again

Revision ID: 883fbe0276f5
Revises: b91baf33e923
Create Date: 2024-11-26 12:02:43.737862

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '883fbe0276f5'
down_revision = 'b91baf33e923'
branch_labels = None
depends_on = None


def upgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('compounds', schema=None) as batch_op:
        batch_op.alter_column('compound_image',
               existing_type=sa.TEXT(),
               type_=sa.String(length=5000),
               existing_nullable=True)

    # ### end Alembic commands ###


def downgrade():
    # ### commands auto generated by Alembic - please adjust! ###
    with op.batch_alter_table('compounds', schema=None) as batch_op:
        batch_op.alter_column('compound_image',
               existing_type=sa.String(length=5000),
               type_=sa.TEXT(),
               existing_nullable=True)

    # ### end Alembic commands ###